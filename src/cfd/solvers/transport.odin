package solvers

import "../../cfd"

transport :: proc(
	mesh: cfd.Mesh,
	system: ^cfd.Linear_System_Builder,
	field: ^cfd.Field,
	mass_flux: cfd.Face_Data,
	diffusivity: f64, density: f64,
	dt: Maybe(f64) = nil,
) -> cfd.Linsolve_Error {
	for _, cell_id in mesh.cells {
		cfd.linear_system_builder_start_cell(system, cell_id)
		defer cfd.linear_system_builder_end_cell(system)
		cfd.add_advection(system, mesh, field^, mass_flux, cell_id)
		if dt, ok := dt.?; ok do cfd.add_time(system, mesh, field^, cell_id, dt / density)
		system.mode = .Sub
		cfd.add_laplacian(system, mesh, field, cell_id, diffusivity)
	}
    solver_mem := cfd.solver_memory_create()
    defer cfd.solver_memory_destroy(&solver_mem)
	ls := cfd.linear_system_builder_assemble(system)
	cfd.linear_system_solve_pgs(ls, field, &solver_mem) or_return
	cfd.linear_system_builder_reset(system)
	return .None
}
