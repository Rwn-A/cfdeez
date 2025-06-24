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

/*    
    The field_gradient field_faces etc family of functions require some care.
    
    If you modify the returned slice the next time that function is called it may return
    the modified data. 
	
	To avoid this, use the detach family of functions, they will clone the data and pass it back to you.
*/
package cfd

import "core:log"
import "core:mem"
import "core:slice"

// It took me a while to decide a global was okay here.
// I decided it was okay because the point of a pool is to reuse memory
// So having a bunch of them doesn't make alot of sense
// Additionally, it used in two unconnected areas, managing fields, and 
// extra buffers for matrix preconditioning. So passing the field pool 
// far enough up the call stack to connect these two didn't seem right.
field_pool: Field_Data_Pool

Boundary_Condition :: struct {
	kind:  enum {
		Dirichlet,
		Neumann,
	},
	value: f64,
}

Boundaries :: map[Boundary_Id]Boundary_Condition

Field :: struct {
	data:                Cell_Data,
	boundary_conditions: Boundaries,
	gradient:            [2]Cell_Data,
	faces:               Face_Data,
	derived_ready:       bit_set[enum {
		Gradient,
		Faces,
	}],
}

Vector_Field :: struct {
	components:    [2]Field,
	divergence:    Cell_Data,
	flux:          Face_Data,
	derived_ready: bit_set[enum {
		Divergence,
		Flux,
	}],
}

Face_Data :: distinct []f64
Cell_Data :: []f64

Primary_Fields :: struct {
	u:        Vector_Field,
	p:        Field,
	passives: []Field,
}

field_set_bnd :: proc(field: ^Field, id: Boundary_Id, bc: Boundary_Condition) {
	field.boundary_conditions[id] = bc
}

field_set_bnd_from_name :: proc(field: ^Field, mesh: Mesh, name: string, bc: Boundary_Condition) {
	field.boundary_conditions[mesh.boundary_names[name]] = bc
}

field_gradient :: proc(mesh: Mesh, field: ^Field) -> [2]Cell_Data {
	if .Gradient not_in field.derived_ready {
		grad := field_pool_vector_cell_data(&field_pool)
		compute_gradient(mesh, field_faces(mesh, field), field.data, grad)
		field.gradient = grad
		field.derived_ready += {.Gradient}
	}
	return field.gradient
}

field_faces :: proc(mesh: Mesh, field: ^Field) -> Face_Data {
	if .Faces not_in field.derived_ready {
		faces_out := field_pool_face_data(&field_pool)
		compute_faces_interpolate(mesh, field^, faces_out)
		field.faces = faces_out
		field.derived_ready += {.Faces}
	}
	return field.faces
}

field_divergence :: proc(mesh: Mesh, vfield: ^Vector_Field) -> Cell_Data {
	if .Divergence not_in vfield.derived_ready {
		div_out := field_pool_cell_data(&field_pool)
		compute_divergence(mesh, field_flux(mesh, vfield), div_out)
		vfield.divergence = div_out
		vfield.derived_ready += {.Divergence}
	}
	return vfield.divergence
}

field_flux :: proc(mesh: Mesh, vfield: ^Vector_Field) -> Face_Data {
	if .Flux not_in vfield.derived_ready {
		flux_out := field_pool_face_data(&field_pool)
		compute_flux(
			mesh,
			{field_faces(mesh, &vfield.components.x), field_faces(mesh, &vfield.components.y)},
			flux_out,
		)
		vfield.flux = flux_out
		vfield.derived_ready += {.Flux}
	}
	return vfield.flux
}

field_dirty :: proc(field: ^Field) {
	if .Gradient in field.derived_ready {
		field_pool_return(&field_pool, field.gradient.x, field.gradient.y)
	}
	if .Faces in field.derived_ready {
		field_pool_return(&field_pool, field.faces)
	}
	field.derived_ready = {}
}

vfield_dirty :: proc(field: ^Vector_Field) {
	field_dirty(&field.components.x)
	field_dirty(&field.components.y)
	if .Divergence in field.derived_ready {
		field_pool_return(&field_pool, field.divergence)
	}
	if .Flux in field.derived_ready {
		field_pool_return(&field_pool, field.flux)
	}
	field.derived_ready = {}
}

field_detach_faces :: proc(mesh: Mesh, field: ^Field) -> Face_Data {
	return field_data_clone(field_faces(mesh, field))
}

field_detach_gradient :: proc(mesh: Mesh, field: ^Field) -> [2]Cell_Data {
	return field_data_clone(field_gradient(mesh, field))
}


vfield_detach_flux :: proc(mesh: Mesh, field: ^Vector_Field) -> Face_Data {
	return field_data_clone(field_flux(mesh, field))
}

vfield_detach_divergence :: proc(mesh: Mesh, field: ^Vector_Field) -> Cell_Data {
	return field_data_clone(field_divergence(mesh, field))
}

cell_data_clone :: proc(cd: Cell_Data) -> Cell_Data {
	field_clone := field_pool_cell_data(&field_pool)
	copy(field_clone, cd)
	return field_clone
}

vector_cell_data_clone :: proc(vcd: [2]Cell_Data) -> [2]Cell_Data {
	return {cell_data_clone(vcd.x), cell_data_clone(vcd.y)}
}

face_data_clone :: proc(faces: Face_Data) -> Face_Data {
	faces_clone := field_pool_face_data(&field_pool)
	copy(faces_clone, faces)
	return faces_clone
}

field_data_clone :: proc {
	cell_data_clone,
	face_data_clone,
	vector_cell_data_clone,
}
