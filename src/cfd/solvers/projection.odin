package solvers

import "../../cfd"
import "core:log"

projection :: proc(
	mesh: cfd.Mesh,
	u: ^cfd.Vector_Field,
	p: ^cfd.Field,
	mass_flux: ^cfd.Face_Data, // pointer pass is not strictly nessecary, just to indicate it's modified.
	systems: [2]^cfd.Linear_System_Builder,
	viscosity, density, dt: f64,
) -> cfd.Linsolve_Error {
	for _, cell_id in mesh.cells {
		cfd.linear_system_builder_start_cell(systems[0], cell_id)
		cfd.linear_system_builder_start_cell(systems[1], cell_id)
		defer cfd.linear_system_builder_end_cell(systems[0])
		defer cfd.linear_system_builder_end_cell(systems[1])
		cfd.add_advection(systems[0], mesh, u.components.x, cfd.field_flux(mesh, u), cell_id)
		cfd.add_advection(systems[1], mesh, u.components.y, cfd.field_flux(mesh, u), cell_id)
		cfd.add_time(systems[0], mesh, u.components.x, cell_id, dt)
		cfd.add_time(systems[1], mesh, u.components.y, cell_id, dt)
		systems[0].mode = .Sub;systems[1].mode = .Sub
		cfd.add_laplacian(systems[0], mesh, &u.components.x, cell_id, viscosity)
		cfd.add_laplacian(systems[1], mesh, &u.components.y, cell_id, viscosity)
	}
    solver_mem := cfd.solver_memory_create()
    defer cfd.solver_memory_destroy(&solver_mem)

	ls_x := cfd.linear_system_builder_assemble(systems[0])
	ls_y := cfd.linear_system_builder_assemble(systems[1])
	cfd.linear_system_solve_pgs(ls_x, &u.components.x, &solver_mem) or_return
	cfd.linear_system_solve_pgs(ls_y, &u.components.y, &solver_mem) or_return
	cfd.linear_system_builder_reset(systems[0])
	cfd.linear_system_builder_reset(systems[1])

	div := cfd.vfield_detach_divergence(mesh, u)
	for &d in div do d *= density / dt

	for cell, cell_id in mesh.cells {
		cfd.linear_system_builder_start_cell(systems[0], cell_id)
		defer cfd.linear_system_builder_end_cell(systems[0])
		cfd.add_laplacian(systems[0], mesh, p, cell_id, 1.0)
	}
	cfd.linear_system_builder_add_source(systems[0], mesh, div)
	ls := cfd.linear_system_builder_assemble(systems[0])
	cfd.linear_system_solve_pcg(ls, p, &solver_mem) or_return
	cfd.linear_system_builder_reset(systems[0])

	grad := cfd.field_gradient(mesh, p)
	for cell, cell_id in mesh.cells {
		u.components.x.data[cell_id] -= (grad.x[cell_id] * dt / density)
		u.components.y.data[cell_id] -= (grad.y[cell_id] * dt / density)
	}
	cfd.vfield_dirty(u)

	// tot := 0.0
	// diverge := cfd.field_divergence(mesh, u)
	// for _, cell_id in mesh.cells {
	//     tot += diverge[cell_id]
	// }
	// log.info(tot / f64(len(mesh.cells)))

	copy(mass_flux^, cfd.field_flux(mesh, u))
	for &d in mass_flux do d *= density

	return .None
}
