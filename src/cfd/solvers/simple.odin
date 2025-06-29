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
import "core:fmt"

// TODO: check for convergence instead of fixed iters.

SIMPLE_ITERS :: 30
MOMENTUM_RELAX :: 1
PRESSURE_RELAX :: 1

SIMPLE_State_Vars :: struct {
	u:       ^cfd.Vector_Field,
	p:       ^cfd.Scalar_Field,
}

SIMPLE :: proc(
	mesh: cfd.Mesh,
	fp: ^cfd.Field_Data_Pool,
	s: SIMPLE_State_Vars,
	systems: [2]^cfd.Linear_System_Builder,
	viscosity: f64,
) -> cfd.Linsolve_Error {
	Ap := cfd.field_pool_cell_data(fp)
	defer cfd.field_pool_return_cell_data(fp, Ap)

	u_old := cfd.field_pool_vector_cell_data(fp)
	defer cfd.field_pool_return_cell_data(fp, ..u_old[:])

	H := cfd.field_pool_vector_cell_data(fp)
	defer cfd.field_pool_return_cell_data(fp, ..H[:])

	p_old := cfd.field_pool_cell_data(fp)
	defer cfd.field_pool_return_cell_data(fp, p_old)

	for iter in 0 ..< SIMPLE_ITERS {
		//setup momentum equation
		for _, cell_id in mesh.cells {
			cfd.linear_system_builder_start_cell(systems[0], cell_id)
			cfd.linear_system_builder_start_cell(systems[1], cell_id)
			defer cfd.linear_system_builder_end_cell(systems[0])
			defer cfd.linear_system_builder_end_cell(systems[1])
			cfd.add_advection(systems[0], mesh, &s.u.components.x, cfd.field_flux(mesh, s.u), cell_id)
			cfd.add_advection(systems[1], mesh, &s.u.components.y, cfd.field_flux(mesh, s.u), cell_id)
			systems[0].mode = .Sub;systems[1].mode = .Sub
			cfd.add_laplacian(systems[0], mesh, &s.u.components.x, cell_id, viscosity)
			cfd.add_laplacian(systems[1], mesh, &s.u.components.y, cell_id, viscosity)
		}
		systems[0].mode = .Sub;systems[1].mode = .Sub
		grad := cfd.field_gradient(mesh, s.p)
		cfd.linear_system_builder_add_source(systems[0], mesh, grad.x)
		cfd.linear_system_builder_add_source(systems[1], mesh, grad.y)
		
		//solve momentum equation and relax
		ls_x := cfd.linear_system_builder_assemble(systems[0])
		ls_y := cfd.linear_system_builder_assemble(systems[1])
		{
			cfd.linear_system_get_diagonal(ls_x, mesh, Ap)

			copy(u_old.x, s.u.components.x.data)
			copy(u_old.y, s.u.components.y.data)

			cfd.field_start_mutation(s.u)
			defer cfd.field_end_mutation(s.u)

			linsolve_memory := cfd.linsolve_memory_from_field_pool(fp)
			defer cfd.linsolve_memory_return_to_pool(&linsolve_memory, fp)

			cfd.linear_system_solve_pgs(ls_x, s.u.components.x.data, &linsolve_memory) or_return
			cfd.linear_system_solve_pgs(ls_y, s.u.components.y.data, &linsolve_memory) or_return

			for _, cell_id in mesh.cells {
				xdat := &s.u.components.x.data[cell_id]
				ydat := &s.u.components.y.data[cell_id]
				xdat^ = MOMENTUM_RELAX * xdat^ + (1.0 - MOMENTUM_RELAX) * u_old.x[cell_id]
				ydat^ = MOMENTUM_RELAX * ydat^ + (1.0 - MOMENTUM_RELAX) * u_old.y[cell_id]
			}
		}

		// build H vector (analogous to velocity but not equal)
		for _, cell_id in mesh.cells {
			H.x[cell_id] = ls_x.b[cell_id]
			H.y[cell_id] = ls_y.b[cell_id]

			for i in ls_x.M.row_indices[cell_id] ..< ls_x.M.row_indices[cell_id + 1] {
				if ls_x.M.columns[i] != cell_id {
					H.x[cell_id] -= ls_x.M.values[i] * s.u.components.x.data[ls_x.M.columns[i]]
				}
			}

			for i in ls_y.M.row_indices[cell_id] ..< ls_y.M.row_indices[cell_id + 1] {
				if ls_y.M.columns[i] != cell_id {
					H.y[cell_id] -= ls_y.M.values[i] * s.u.components.y.data[ls_y.M.columns[i]]
				}
			}

			H.x[cell_id] /= Ap[cell_id]
			H.y[cell_id] /= Ap[cell_id]
		}

		// TODO: these BC's aren't completely accurate at outflow's because of non-zero pressure gradient
		// TODO: helper func for this sort of thing
		temp_H := cfd.Vector_Field {
			components     = {
				{data = H.x, bnds = s.u.components.x.bnds, pool = fp, _derived_ready = {}},
				{data = H.y, bnds = s.u.components.y.bnds, pool = fp, _derived_ready = {}},
			},
			_derived_ready = {},
			pool           = fp,
		}
		div := cfd.field_divergence(mesh, &temp_H)
		cfd.linear_system_builder_reset(systems[0])
		cfd.linear_system_builder_reset(systems[1])

		copy(p_old, s.p.data)
		//clear pressure, its a correction here.
		cfd.scalar_field_multiply_scalar(s.p, 0)

		//setup pressure equation
		for cell, cell_id in mesh.cells {
			cfd.linear_system_builder_start_cell(systems[0], cell_id)
			defer cfd.linear_system_builder_end_cell(systems[0])
			cfd.add_laplacian(systems[0], mesh, s.p, cell_id, 1.0 / Ap[cell_id])
		}
		cfd.linear_system_builder_add_source(systems[0], mesh, div)

		//solve pressure equation
		ls_p := cfd.linear_system_builder_assemble(systems[0])
		{
			cfd.field_start_mutation(s.p)
			defer cfd.field_end_mutation(s.p)

			linsolve_memory := cfd.linsolve_memory_from_field_pool(fp)
			defer cfd.linsolve_memory_return_to_pool(&linsolve_memory, fp)

			cfd.linear_system_solve_pcg(ls_p, s.p.data, &linsolve_memory) or_return
			cfd.linear_system_builder_reset(systems[0])
		}

		// correct velocities & update pressure
		{
			cfd.field_start_mutation(s.u)
			defer cfd.field_end_mutation(s.u)

			grad := cfd.field_gradient(mesh, s.p)
			for _, cell_id in mesh.cells {
				d := 1.0 / Ap[cell_id]
				s.u.components.x.data[cell_id] -= d * grad.x[cell_id]
				s.u.components.y.data[cell_id] -= d * grad.y[cell_id]
			}

			cfd.field_start_mutation(s.p)
			defer cfd.field_end_mutation(s.p)

			for _, cell_id in mesh.cells {
				s.p.data[cell_id] = p_old[cell_id] + PRESSURE_RELAX * s.p.data[cell_id]
			}
		}

	}
	return .None
}
