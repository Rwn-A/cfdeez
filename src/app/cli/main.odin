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

run :: proc(config_path: string) -> bool {

	log.info("Loading config...")
	c: app.Case
	app.load_case(&c, config_path) or_return
	defer vmem.arena_destroy(&c.arena)

	log.info("Preparing Simulation...")

	if !os.exists(c.output.directory) {
		if err := os.make_directory(c.output.directory); err != nil {
			log.errorf("Unable to create output directory %s", c.output.directory)
			return false
		}
	}

	arena_alloc := vmem.arena_allocator(&c.arena)

	passive_system := cfd.linear_system_builder_new(c.mesh, arena_alloc)
	u_system_x := cfd.linear_system_builder_new(c.mesh, arena_alloc) //pressure can reuse one of these.
	u_system_y := cfd.linear_system_builder_new(c.mesh, arena_alloc)

	_, is_transient := c.timestep.?

	write_pvd := .VTK in c.output.format && is_transient
	pvd_writer: cfd_io.VTK_Writer
	if write_pvd {
		pvd_writer = cfd_io.vtk_start_pvd(c.output.directory, c.name) or_return
	}
	defer if write_pvd {
		cfd_io.vtk_end_pvd(&pvd_writer)
	}

	output_scalar_fields := slice.concatenate([][]cfd.Scalar_Field{[]cfd.Scalar_Field{c.fields.p}, c.fields.passives})
	output_scalar_names := slice.concatenate([][]string{[]string{"pressure"}, c.passive_names})
	defer delete(output_scalar_fields)
	defer delete(output_scalar_names)

	output :: proc(c: app.Case, sfs: []cfd.Scalar_Field, sf_names: []string, step: int) -> (vtk_name: string, ok: bool) {
		dir := c.output.directory
		for opt in c.output.format {
			switch opt {
			case .VTK:
				vtk_name = cfd_io.vtk_write_vtu(
					c.mesh,
					{c.fields.u},
					sfs,
					{"velocity"},
					sf_names,
					dir,
					c.name,
					step,
				) or_return
			case .CSV:
				_ = cfd_io.write_csv(c.mesh, {c.fields.u}, sfs, {"velocity"}, sf_names, dir, c.name, step) or_return
			}
		}
		return vtk_name, true
	}

	if vtk_name, has := output(c, output_scalar_fields, output_scalar_names, 0); has && write_pvd {
		cfd_io.vtk_write_pvd_entry(&pvd_writer, vtk_name, 0)
	}

	update_mass_flux :: proc(mesh: cfd.Mesh, u: ^cfd.Vector_Field, flux: cfd.Face_Data, density: f64) {
		volume_flux := cfd.field_flux(mesh, u)
		for &d, i in flux {d = volume_flux[i] * density}
	}

	mass_flux := cfd.field_pool_face_data(&c.fp)
	update_mass_flux(c.mesh, &c.fields.u, mass_flux, c.fluid.density)
	log.info("Simulation running...")

	for i in 1 ..= c.steps {
		if .IncFlow in c.physics && is_transient {
			err := solvers.pressure_projection(
				c.mesh,
				&c.fp,
				{u = &c.fields.u, p = &c.fields.p},
				{&u_system_x, &u_system_y},
				c.fluid.viscosity,
				c.timestep.?,
			)
			if err != .None {
				log.info("Simulation failed.");return false
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
				log.info("Simulation failed.");return false
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
					log.info("Simulation failed.");return false
				}
			}
		}
		if i % c.output.frequency == 0 {
			if vtk_name, ok := output(c, output_scalar_fields, output_scalar_names, i); ok && write_pvd {
				cfd_io.vtk_write_pvd_entry(&pvd_writer, vtk_name, f64(i) * c.timestep.?)
			}
		}
	}
	log.infof("Simulation complete! Output written to %s", c.output.directory)

	return true
}
