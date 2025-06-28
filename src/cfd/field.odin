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
package cfd

import "core:mem"
import "core:mem/virtual"
import "core:slice"

VF_COMPONENTS :: 2

Boundary_Condition :: struct {
	kind:  enum {
		Dirichlet,
		Neumann,
	},
	value: f64,
}

Boundaries :: map[BoundaryId]Boundary_Condition

// distinct mainly to avoid accidently passing face data into things
// if field data was also distinct it would be more of a pain because the spla code assumes []f64.
Face_Data :: distinct []f64

Scalar_Field :: struct {
	bnds:           Boundaries,
	data:           []f64,
	pool:           ^Field_Data_Pool,
	_gradient:      [VF_COMPONENTS][]f64,
	_faces:         Face_Data,
	_derived_ready: bit_set[enum {
		Gradient,
		Faces,
	}],
	_modifiable:    bool,
}

Vector_Field :: struct {
	components:     [VF_COMPONENTS]Scalar_Field,
	pool:           ^Field_Data_Pool,
	_divergence:    []f64,
	_flux:          Face_Data,
	_derived_ready: bit_set[enum {
		Divergence,
		Flux,
	}],
}

scalar_field_new :: proc(pool: ^Field_Data_Pool, bnd_allocator := context.allocator) -> Scalar_Field {
	return {data = field_pool_cell_data(pool), bnds = make(Boundaries, bnd_allocator), pool = pool}
}

vector_field_new :: proc(pool: ^Field_Data_Pool, bnd_allocator := context.allocator) -> (vf: Vector_Field) {
	vf.components.x = scalar_field_new(pool, bnd_allocator)
	vf.components.y = scalar_field_new(pool, bnd_allocator)
	vf.pool = pool
	return vf
}

vector_field_at :: proc(vf: Vector_Field, idx: CellId) -> Vector {
	return {vf.components.x.data[idx], vf.components.y.data[idx]}
}

/*
	Note for future me: Whats the point of all this?
	Basically, I wan't to do lazy computation for derived fields, after a mutation they need to be marked dirty
	Instead of relying on the user to know when things should be marked dirty we protect the memory.
	In release mode, we take off the protections for performance. This pattern is also pretty natural.

	set up some equations by reading fields and derived fields
	solve equation
	start mutation
	modify field based on result
	end mutation
	read out derived fields again for whatever else.
*/

// -- Mutation control --

@(private)
scalar_field_start_mutation :: proc(sf: ^Scalar_Field) {
	sf._modifiable = true
	if .Faces in sf._derived_ready { field_pool_return_face_data(sf.pool, sf._faces); sf._faces = {} }
	if .Gradient in sf._derived_ready { field_pool_return_cell_data(sf.pool, sf._gradient.x) }
	if .Gradient in sf._derived_ready { field_pool_return_cell_data(sf.pool, sf._gradient.y) }
	sf._gradient = {}
	sf._derived_ready = {}
}

@(private)
scalar_field_end_mutation :: proc(sf: ^Scalar_Field) {
	sf._modifiable = false
}

@(private)
vector_field_start_mutation :: proc(vf: ^Vector_Field) {
	scalar_field_start_mutation(&vf.components.x)
	scalar_field_start_mutation(&vf.components.y)
	if .Flux in vf._derived_ready { field_pool_return_face_data(vf.pool, vf._flux); vf._flux = {} }
	if .Divergence in vf._derived_ready { field_pool_return_cell_data(vf.pool, vf._divergence); vf._divergence = {} }
	vf._derived_ready = {}
}

@(private)
vector_field_end_mutation :: proc(vf: ^Vector_Field) {
	scalar_field_end_mutation(&vf.components.x)
	scalar_field_end_mutation(&vf.components.y)
}

// deletes all derived fields, disables recomputation of derived fields.
field_start_mutation :: proc {
	scalar_field_start_mutation,
	vector_field_start_mutation,
}

// re enables access to derived fields, safe to call even if no modification were started if that need arises.
field_end_mutation :: proc {
	vector_field_end_mutation,
	scalar_field_end_mutation,
}

// -- Derived fields --

// Note: not memory protecting these because modifying them might make sense to avoid a clone.

field_gradient :: proc(mesh: Mesh, field: ^Scalar_Field) -> [VF_COMPONENTS][]f64 {
	assert(
		!field._modifiable,
		"Cannot access field gradient when mutation mode is active. Call field_end_mutation first.",
	)
	if .Gradient not_in field._derived_ready {
		grad := field_pool_vector_cell_data(field.pool)
		compute_gradient(mesh, field_faces(mesh, field), field.data, grad)
		field._gradient = grad
		field._derived_ready += {.Gradient}
	}
	return field._gradient
}

field_faces :: proc(mesh: Mesh, field: ^Scalar_Field) -> Face_Data {
	assert(
		!field._modifiable,
		"Cannot access field faces when mutation mode is active. Call field_end_mutation first.",
	)
	if .Faces not_in field._derived_ready {
		faces_out := field_pool_face_data(field.pool)
		compute_faces_interpolate(mesh, field^, faces_out)
		field._faces = faces_out
		field._derived_ready += {.Faces}
	}
	return field._faces
}

field_divergence :: proc(mesh: Mesh, field: ^Vector_Field) -> []f64 {
	assert(
		!field.components.x._modifiable && !field.components.y._modifiable,
		"Cannot access field flux when mutation mode is active. Call field_end_mutation first.",
	)
	if .Divergence not_in field._derived_ready {
		div_out := field_pool_cell_data(field.pool)
		compute_divergence(mesh, field_flux(mesh, field), div_out)
		field._divergence = div_out
		field._derived_ready += {.Divergence}
	}
	return field._divergence
}

field_flux :: proc(mesh: Mesh, field: ^Vector_Field) -> Face_Data {
	assert(
		!field.components.x._modifiable && !field.components.y._modifiable,
		"Cannot access field divergence when mutation mode is active. Call field_end_mutation first.",
	)
	if .Flux not_in field._derived_ready {
		flux_out := field_pool_face_data(field.pool)
		compute_flux(
			mesh,
			[VF_COMPONENTS]Face_Data{field_faces(mesh, &field.components.x), field_faces(mesh, &field.components.y)},
			flux_out,
		)
		field._flux = flux_out
		field._derived_ready += {.Flux}
	}
	return field._flux
}


// -- Some helpers to avoid constant start, stop modifiers -- 

scalar_field_multiply_scalar :: proc(sf: ^Scalar_Field, scalar: f64) {
	field_start_mutation(sf)
	for &d in sf.data {d *= scalar}
	field_end_mutation(sf)
}

vector_field_multiply_scalar :: proc(vf: ^Vector_Field, scalar: f64) {
	for &component in vf.components {scalar_field_multiply_scalar(&component, scalar)}
}

field_multiply_scalar :: proc {
	scalar_field_multiply_scalar,
	vector_field_multiply_scalar,
}

// -- Field pool --

Field_Data_Pool :: struct {
	backing:              mem.Allocator,
	num_faces, num_cells: int,
	free_cell_data:       [dynamic][]f64,
	free_face_data:       [dynamic]Face_Data,
}

field_pool_init :: proc(pool: ^Field_Data_Pool, mesh: Mesh, backing := context.allocator) {
	DEFAULT_ALLOCATED_FIELDS :: 8

	pool.backing = backing
	context.allocator = backing

	pool.num_cells = len(mesh.cells)
	pool.num_faces = len(mesh.faces)

	pool.free_cell_data = make([dynamic][]f64, 0, DEFAULT_ALLOCATED_FIELDS)
	pool.free_face_data = make([dynamic]Face_Data, 0, DEFAULT_ALLOCATED_FIELDS)

	for _ in 0 ..< DEFAULT_ALLOCATED_FIELDS {
		append(&pool.free_cell_data, make([]f64, pool.num_cells))
		append(&pool.free_face_data, make(Face_Data, pool.num_faces))
	}
}

// This will only destroy fields that have been returned.
// If all fields are not returned before this call, they will leak unless the backing allocator is freed.
field_pool_destroy :: proc(pool: ^Field_Data_Pool) {
	for f in pool.free_cell_data {delete(f, pool.backing)}
	for f in pool.free_face_data {delete(f, pool.backing)}
	delete(pool.free_cell_data)
	delete(pool.free_face_data)
}

field_pool_cell_data :: proc(pool: ^Field_Data_Pool) -> []f64 {
	data := pop_safe(&pool.free_cell_data) or_else make([]f64, pool.num_cells, pool.backing)
	return data
}

field_pool_face_data :: proc(pool: ^Field_Data_Pool) -> Face_Data {
	data := pop_safe(&pool.free_face_data) or_else make(Face_Data, pool.num_faces, pool.backing)
	return data
}

field_pool_vector_cell_data :: proc(pool: ^Field_Data_Pool) -> [VF_COMPONENTS][]f64 {
	return {field_pool_cell_data(pool), field_pool_cell_data(pool)}
}

field_pool_return_cell_data :: proc(pool: ^Field_Data_Pool, items: ..[]f64) {
	for cell_data in items {
		when ODIN_DEBUG {
			slice.fill(cell_data, 0)
		}
		assert(len(cell_data) == pool.num_cells)
		append(&pool.free_cell_data, cell_data)
	}
}

field_pool_return_face_data :: proc(pool: ^Field_Data_Pool, items: ..Face_Data) {
	for face_data in items {
		when ODIN_DEBUG {slice.fill(face_data, 0)}
		assert(len(face_data) == pool.num_faces)
		append(&pool.free_face_data, face_data)
	}
}

@(private)
cell_data_clone :: proc(pool: ^Field_Data_Pool, data: []f64) -> []f64 {
	field_clone := field_pool_cell_data(pool)
	copy(field_clone, data)
	return field_clone
}

@(private)
vector_cell_data_clone :: proc(pool: ^Field_Data_Pool, data: [VF_COMPONENTS][]f64) -> [VF_COMPONENTS][]f64 {
	return {cell_data_clone(pool, data.x), cell_data_clone(pool, data.y)}
}

@(private)
face_data_clone :: proc(pool: ^Field_Data_Pool, data: Face_Data) -> Face_Data {
	faces_clone := field_pool_face_data(pool)
	copy(faces_clone, data)
	return faces_clone
}

field_data_clone :: proc {
	cell_data_clone,
	face_data_clone,
	vector_cell_data_clone,
}
