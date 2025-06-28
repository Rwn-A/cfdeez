/*
 Copyright (C) 2025 Rowan Apps, Tor Rabien

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
package main

import "core:log"
import vmem "core:mem/virtual"
import "core:os"
import "core:slice"
import "core:thread"
import "core:sync"

import "../../app"
import "../../cfd"
import cfd_io "../../cfd/io"
import solvers "../../cfd/solvers"

main :: proc() {
	context.logger = log.create_console_logger()
	defer log.destroy_console_logger(context.logger)

	if len(os.args) < 2 {
		log.error("Expected a path to a config file.")
		os.exit(1)
	}

	if ok := run(os.args[1]); !ok {
		os.exit(1)
	}
}

// output is on a seperate thread. This only makes a difference for transient sims with lots file writes
// I'm not experienced with this sort of thing, from what i can gather this should be okay?. 
// If the output is really slow then maybe a I'll add a queue or something

Thread_State :: struct {
	buffers: ^Output_Buffers,
	step: ^int,
	done: ^bool,
	ready: ^sync.Cond,
}

Output_Buffers :: struct {
	out_u: cfd_io.Out_Vector_Field,
	out_scalars: []cfd_io.Out_Scalar_Field, // index 0 is pressure.

	mutex: sync.Mutex,
}

run :: proc(config_path: string) -> bool {
	log.info("Loading config...")

	c: app.Case
	app.load_case(&c, config_path) or_return
	defer vmem.arena_destroy(&c.arena)
	
	log.info("Preparing Simulation...")

	arena_alloc := vmem.arena_allocator(&c.arena)
	passive_system := cfd.linear_system_builder_new(c.mesh, arena_alloc)
	u_system_x := cfd.linear_system_builder_new(c.mesh, arena_alloc) //pressure can reuse one of these.
	u_system_y := cfd.linear_system_builder_new(c.mesh, arena_alloc)

	//prepare buffers for the output thread.
	output_buffers := Output_Buffers{
		out_u = {cfd.field_pool_vector_cell_data(&c.fp), "velocity"},
	}
	out_scalars := make([dynamic]cfd_io.Out_Scalar_Field)
	defer delete(out_scalars)

	append(&out_scalars, cfd_io.Out_Scalar_Field{cfd.field_pool_cell_data(&c.fp), "pressure"})
	for name in c.passive_names{
		append(&out_scalars, cfd_io.Out_Scalar_Field{cfd.field_pool_cell_data(&c.fp), name})
	}
	output_buffers.out_scalars = out_scalars[:]

	done: bool
	output_ready := sync.Cond{}
	step := 0

	output_thread_state := Thread_State{
		buffers = &output_buffers,
		done = &done,
		ready = &output_ready,
		step = &step
	}

	output_thread := thread.create_and_start_with_poly_data2(&output_thread_state, &c, write_output)

	// output initial conditions 
	sync.lock(&output_buffers.mutex)
    write_output_buffers(output_buffers, c)
    sync.cond_signal(&output_ready)
    sync.unlock(&output_buffers.mutex)

	dt, is_transient := c.timestep.?

	update_mass_flux :: proc(mesh: cfd.Mesh, u: ^cfd.Vector_Field, flux: cfd.Face_Data, density: f64) {
		volume_flux := cfd.field_flux(mesh, u)
		for &d, i in flux {d = volume_flux[i] * density}
	}

	mass_flux := cfd.field_pool_face_data(&c.fp)
	update_mass_flux(c.mesh, &c.fields.u, mass_flux, c.fluid.density)
	log.info("Simulation running...")

	for i in 1 ..= c.steps {
		step = i
		if .IncFlow in c.physics && is_transient {
			err := solvers.pressure_projection(
				c.mesh,
				&c.fp,
				{u = &c.fields.u, p = &c.fields.p},
				{&u_system_x, &u_system_y},
				c.fluid.viscosity,
				dt,
			)
			if err != .None {
				log.info("Simulation failed."); done = true; return false
			}
		}
		if .IncFlow in c.physics && !is_transient {
			err := solvers.SIMPLE(
				c.mesh,
				&c.fp,
				{u = &c.fields.u, p = &c.fields.p},
				{&u_system_x, &u_system_y},
				c.fluid.viscosity,
			)
			if err != .None {
				log.info("Simulation failed."); done = true; return false
			}
		}
		update_mass_flux(c.mesh, &c.fields.u, mass_flux, c.fluid.density)
		if .Transport in c.physics {
			for &passive, idx in c.fields.passives {
				err := solvers.transport(
					c.mesh,
					&c.fp,
					{scalar = &passive, density = c.fluid.density},
					&passive_system,
					mass_flux,
					c.passive_diffusivities[idx],
					c.timestep,
				)
				if err != .None {
					log.info("Simulation failed."); done = true; return false
				}
			}
		}
		if i % c.output.frequency == 0 {
			sync.lock(&output_buffers.mutex)
            write_output_buffers(output_buffers, c)
            sync.cond_signal(&output_ready)
            sync.unlock(&output_buffers.mutex)
		}
	}
	log.infof("Simulation complete! Output written to %s", c.output.directory)


	sync.lock(&output_buffers.mutex)
    done = true
    sync.cond_signal(&output_ready)
    sync.unlock(&output_buffers.mutex)
    
    thread.join(output_thread)

	return true
}

write_output_buffers :: proc(buffers: Output_Buffers, c: app.Case) {
	copy(buffers.out_u.data.x, c.fields.u.components.x.data)
	copy(buffers.out_u.data.y, c.fields.u.components.y.data)
	copy(buffers.out_scalars[0].data, c.fields.p.data)
	for i in 1..<len(buffers.out_scalars) {
		copy(buffers.out_scalars[i].data, c.fields.passives[i - 1].data)
	}
}

//Note: I know we pass in the case, but I do not access the fields from that, big ol race condition.
write_output :: proc(state: ^Thread_State, 	c: ^app.Case) {
	pvd_writer: cfd_io.VTK_Writer

	dt, is_transient := c.timestep.?
	if is_transient {
		pvd_writer, _ = cfd_io.vtk_start_pvd(c.output.directory, c.name)
	}
	defer if is_transient {
		cfd_io.vtk_end_pvd(&pvd_writer)
	}

	for !state.done^{
		sync.lock(&state.buffers.mutex)
		sync.cond_wait(state.ready, &state.buffers.mutex)

		if !state.done^ {
			if .VTK in c.output.format{
				path := cfd_io.output_path(c.output.directory, c.name, "vtu", state.step^)
				if ok := cfd_io.vtk_write_vtu(c.mesh, {state.buffers.out_u}, state.buffers.out_scalars, path); ok {
					cfd_io.vtk_write_pvd_entry(&pvd_writer, path, dt * f64(state.step^))
				}
			}
			if .CSV in c.output.format{
				path := cfd_io.output_path(c.output.directory, c.name, "csv", state.step^)
				cfd_io.write_csv(c.mesh, {state.buffers.out_u}, state.buffers.out_scalars, path) 
			}
		}
		sync.unlock(&state.buffers.mutex)
	}	
}
