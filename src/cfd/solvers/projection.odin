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
package solvers

import "../../cfd"

Projection_State_Vars :: struct {
	u:       ^cfd.Vector_Field,
	p:       ^cfd.Scalar_Field,
}


pressure_projection :: proc(
	mesh: cfd.Mesh,
	fp: ^cfd.Field_Data_Pool,
	s: Projection_State_Vars,
	systems: [2]^cfd.Linear_System_Builder,
	viscosity: f64,
	dt: f64,
) -> cfd.Linsolve_Error {
	//build momentum equation
	for _, cell_id in mesh.cells {
		cfd.linear_system_builder_start_cell(systems[0], cell_id)
		cfd.linear_system_builder_start_cell(systems[1], cell_id)
		defer cfd.linear_system_builder_end_cell(systems[0])
		defer cfd.linear_system_builder_end_cell(systems[1])
		cfd.add_advection(systems[0], mesh, &s.u.components.x, cfd.field_flux(mesh, s.u), cell_id)
		cfd.add_advection(systems[1], mesh, &s.u.components.y, cfd.field_flux(mesh, s.u), cell_id)
		cfd.add_time(systems[0], mesh, s.u.components.x, cell_id, dt)
		cfd.add_time(systems[1], mesh, s.u.components.y, cell_id, dt)
		systems[0].mode = .Sub;systems[1].mode = .Sub
		cfd.add_laplacian(systems[0], mesh, &s.u.components.x, cell_id, viscosity)
		cfd.add_laplacian(systems[1], mesh, &s.u.components.y, cell_id, viscosity)
	}
	//solve momentum equation
	{
		//we want to do this first so the memory from the derived fields is available to be used as solver memory
		cfd.field_start_mutation(s.u)
		defer cfd.field_end_mutation(s.u)

		linsolve_memory := cfd.linsolve_memory_from_field_pool(fp)
		defer cfd.linsolve_memory_return_to_pool(&linsolve_memory, fp)

		ls_x := cfd.linear_system_builder_assemble(systems[0])
		ls_y := cfd.linear_system_builder_assemble(systems[1])
		cfd.linear_system_solve_pgs(ls_x, s.u.components.x.data, &linsolve_memory) or_return
		cfd.linear_system_solve_pgs(ls_y, s.u.components.y.data, &linsolve_memory) or_return

		cfd.linear_system_builder_reset(systems[0])
		cfd.linear_system_builder_reset(systems[1])
	}
	//build pressure equation

	//modifying divergence is not really a good idea,
	//but it is not used before u is marked modified again, so in this case it's fine.
	div := cfd.field_divergence(mesh, s.u)
	for &d in div {d *= 1.0 / dt}

	for cell, cell_id in mesh.cells {
		cfd.linear_system_builder_start_cell(systems[0], cell_id)
		defer cfd.linear_system_builder_end_cell(systems[0])
		cfd.add_laplacian(systems[0], mesh, s.p, cell_id, 1.0)
	}
	cfd.linear_system_builder_add_source(systems[0], mesh, div)

	//solve pressure equation
	{
		cfd.field_start_mutation(s.p)
		defer cfd.field_end_mutation(s.p)

		linsolve_memory := cfd.linsolve_memory_from_field_pool(fp)
		defer cfd.linsolve_memory_return_to_pool(&linsolve_memory, fp)

		ls := cfd.linear_system_builder_assemble(systems[0])
		cfd.linear_system_solve_pcg(ls, s.p.data, &linsolve_memory) or_return
		cfd.linear_system_builder_reset(systems[0])
	}

	//correct velocities 
	{
		cfd.field_start_mutation(s.u)
		defer cfd.field_end_mutation(s.u)

		grad := cfd.field_gradient(mesh, s.p)
		for cell, cell_id in mesh.cells {
			s.u.components.x.data[cell_id] -= (grad.x[cell_id] * dt )
			s.u.components.y.data[cell_id] -= (grad.y[cell_id] * dt )
		}
	}

	return .None
}
