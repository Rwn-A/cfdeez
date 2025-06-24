package main

import "core:log"
import vmem "core:mem/virtual"
import "core:os"

import "../../cfd"
import cfd_io "../../cfd/io"
import solvers "../../cfd/solvers"
import shared "../../cmd"

main :: proc() {
	context.logger = log.create_console_logger()
	defer log.destroy_console_logger(context.logger)

	if len(os.args) < 2 {
		log.error("Expected a path to a config file.")
		os.exit(1)
	}

	arena: vmem.Arena
	c, case_ok := shared.load_case(os.args[1], &arena)
	if !case_ok do os.exit(1)
	defer vmem.arena_destroy(c.arena)

	if !os.exists(c.output.directory) {
		if err := os.make_directory(c.output.directory); err != nil {
			log.errorf("Unable to create output directory %s", c.output.directory)
			return
		}
	}

	// Note: mass flux makes sense to be outside any particular solver
	// for scalar transport, it only needs that not velocity
	// and for SIMPLE / PISO, mass flux is a seperate quantity not recomputed from U each time.
	mass_flux := cfd.vfield_detach_flux(c.mesh, &c.fields.u)

	for &d in mass_flux do d *= c.fluid.density

	passive_system := cfd.linear_system_builder_new(c.mesh, vmem.arena_allocator(c.arena))

	u_system_x := cfd.linear_system_builder_new(c.mesh, vmem.arena_allocator(c.arena))
	u_system_y := cfd.linear_system_builder_new(c.mesh, vmem.arena_allocator(c.arena))

	pvd_w: cfd_io.VTU_Writer
	if c.output.pvd {
		output_path := cfd_io.output_path(c.output.directory, c.name, "pvd", 0)
		if fd, err := os.open(output_path, os.O_CREATE | os.O_TRUNC | os.O_RDWR, 0o777); err == nil {
			pvd_w = cfd_io.write_pvd_header(fd)
		} else {
			log.errorf("Could not create output PVD file %s", output_path)
			return
		}
	}
	defer if c.output.pvd {
		cfd_io.write_pvd_end(&pvd_w)
		os.close(pvd_w.fd)
	}

	for output in c.output.format do output(c.mesh, c.fields, c.passive_names, c.output.directory, c.name, 0)
	if c.output.pvd {
		cfd_io.write_pvd_entry(&pvd_w, cfd_io.output_path(".", c.name, "vtu", 0), 0)
	}

	log.info("Simulation Running...")

	for i in 1 ..= c.steps {
		if .IncFlow in c.physics {
			if _, has := c.timestep.?; !has {
				err := solvers.simple(
					c.mesh,
					&c.fields.u,
					&c.fields.p,
					&mass_flux,
					{&u_system_x, &u_system_y},
					c.fluid.viscosity,
					c.fluid.density,
				)
				if err != nil {log.errorf("simulation failed. %v", err);return}
			}else{
				err := solvers.projection(
				c.mesh,
				&c.fields.u,
				&c.fields.p,
				&mass_flux,
				{&u_system_x, &u_system_y},
				c.fluid.viscosity,
				c.fluid.density,
				c.timestep.?,
				)
				if err != nil {log.errorf("simulation failed. %v", err);return}
			}
			
		}

		if .Transport in c.physics {
			for &passive, i in c.fields.passives {
				err := solvers.transport(
					c.mesh,
					&passive_system,
					&passive,
					mass_flux,
					c.passive_diffusivity[i],
					c.fluid.density,
					c.timestep,
				)
				if err != nil {log.errorf("simulation failed. %v", err);return}
			}
		}

		if i % c.output.frequency == 0 {
			for output in c.output.format do output(c.mesh, c.fields, c.passive_names, c.output.directory, c.name, i)
			if !c.output.pvd do continue
			cfd_io.write_pvd_entry(&pvd_w, cfd_io.output_path(".", c.name, "vtu", i), f64(i) * c.timestep.?)
		}
	}

	log.infof("Simulation complete! Output written to %s", c.output.directory)
}
