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

Transport_State_Vars :: struct {
	scalar:  ^cfd.Scalar_Field,
	density: f64,
}

transport :: proc(
	mesh: cfd.Mesh,
	fp: ^cfd.Field_Data_Pool,
	s: Transport_State_Vars,
	system: ^cfd.Linear_System_Builder,
	mass_flux: cfd.Face_Data,
	diffusivity: f64,
	dt: Maybe(f64) = nil,
) -> cfd.Linsolve_Error {
	for _, cell_id in mesh.cells {
		cfd.linear_system_builder_start_cell(system, cell_id)
		defer cfd.linear_system_builder_end_cell(system)
		cfd.add_advection(system, mesh, s.scalar^, mass_flux, cell_id)
		if dt, ok := dt.?; ok do cfd.add_time(system, mesh, s.scalar^, cell_id, dt / s.density)
		system.mode = .Sub
		cfd.add_laplacian(system, mesh, s.scalar, cell_id, diffusivity)
	}
	{
		cfd.field_start_mutation(s.scalar)
		defer cfd.field_end_mutation(s.scalar)

		linsolve_memory := cfd.linsolve_memory_from_field_pool(fp)
		defer cfd.linsolve_memory_return_to_pool(&linsolve_memory, fp)

		ls := cfd.linear_system_builder_assemble(system)

		cfd.linear_system_solve_pgs(ls, s.scalar.data, &linsolve_memory) or_return
		cfd.linear_system_builder_reset(system)
	}
	
	return .None
}
