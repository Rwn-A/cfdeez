package solvers

import "../../cfd"

transport :: proc(
    mesh: cfd.Mesh,
    system: ^cfd.Linear_System_Builder,
    field: ^cfd.Field,
    mass_flux: cfd.Face_View,
    diffusivity: f64,
    dt: Maybe(f64)
) -> cfd.Linsolve_Error {
    for _, cell_id in mesh.cells {
        cfd.linear_system_builder_start_cell(system, cell_id)
        defer cfd.linear_system_builder_end_cell(system)
        cfd.add_advection(system, mesh, field^, mass_flux, cell_id)
        if dt, ok := dt.?; ok do cfd.add_time(system, mesh, field^, cell_id, dt)
        system.mode = .Sub
        cfd.add_laplacian(system, mesh, field, cell_id, diffusivity)
    }
    ls := cfd.linear_system_builder_assemble(system)
    cfd.linear_system_solve_gs(ls, field) or_return
    cfd.linear_system_builder_reset(system)
    return .None
}