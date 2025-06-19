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
import "core:slice"
import "core:log"

Boundary_Condition :: struct { kind: enum {Dirichlet, Neumann}, value: f64 }
Boundaries :: map[Boundary_Id]Boundary_Condition

Field :: struct {
    data: Field_View,
    boundary_conditions: Boundaries,
    gradient: Vector_Field_View,
    faces: Face_View,
    derived_ready: bit_set[enum{Gradient, Faces}],
    pool: ^Field_View_Pool,
}

Vector_Field :: struct {
    components: [2]Field,
    divergence: Field_View,
    flux: Face_View,
    derived_ready: bit_set[enum{Divergence, Flux}],
    pool: ^Field_View_Pool,
}

Face_View :: distinct []f64
Field_View :: distinct []f64
Vector_Field_View :: distinct [2]Field_View

Primary_Fields :: struct {
    u: Vector_Field,
    p: Field,
    passives: []Field,
}

field_set_bnd :: proc(field: ^Field, id: Boundary_Id, bc: Boundary_Condition) {
    field.boundary_conditions[id] = bc
}

field_set_bnd_from_name :: proc(field: ^Field, mesh: Mesh, name: string, bc: Boundary_Condition) {
    field.boundary_conditions[mesh.boundary_names[name]] = bc
}

field_gradient :: proc(mesh: Mesh, field: ^Field) -> Vector_Field_View {
    if .Gradient not_in field.derived_ready {
        grad := field_pool_acquire_vfield(field.pool)
        compute_gradient(mesh, field_faces(mesh, field), grad)
        field.gradient = grad
        field.derived_ready += {.Gradient}
    }
    return field.gradient
}

field_faces :: proc(mesh: Mesh, field: ^Field) -> Face_View {
    if .Faces not_in field.derived_ready {
        faces_out := field_pool_acquire_faces(field.pool)
        compute_faces_interpolate(mesh, field^, faces_out)
        field.faces = faces_out
        field.derived_ready += {.Faces}
    }
    return field.faces
}

field_divergence :: proc(mesh: Mesh, vfield: ^Vector_Field) -> Field_View {
    if .Divergence not_in vfield.derived_ready {
        div_out := field_pool_acquire_field(vfield.pool)
        compute_divergence(mesh, field_flux(mesh, vfield), div_out)
        vfield.divergence = div_out
        vfield.derived_ready += {.Divergence}
    }
    return vfield.divergence
}

field_flux :: proc(mesh: Mesh, vfield: ^Vector_Field) -> Face_View {
    if .Flux not_in vfield.derived_ready {
        flux_out := field_pool_acquire_faces(vfield.pool)
        compute_flux(mesh, {field_faces(mesh, &vfield.components.x), field_faces(mesh, &vfield.components.y)}, flux_out)
        vfield.flux = flux_out
        vfield.derived_ready += {.Flux}
    }
    return vfield.flux
}

field_dirty :: proc(field: ^Field) {
    if .Gradient in field.derived_ready {
        field_pool_return(field.pool, field.gradient.x, field.gradient.y)
    }
    if .Faces in field.derived_ready {
        field_pool_return(field.pool, field.faces)
    }
    field.derived_ready = {}
}

vfield_dirty :: proc(field: ^Vector_Field) {
    field_dirty(&field.components.x)
    field_dirty(&field.components.y)
    if .Divergence in field.derived_ready {
        field_pool_return(field.pool, field.divergence)
    }
    if .Flux in field.derived_ready {
        field_pool_return(field.pool, field.flux)
    }
    field.derived_ready = {}
}

Field_View_Pool :: struct {
    num_faces, num_cells: int,
    free_fields: [dynamic]Field_View,
    free_faces: [dynamic]Face_View,
    backing_allocator: mem.Allocator,
}

field_view_clone :: proc(field: Field_View, pool: ^Field_View_Pool) -> Field_View {
    field_clone := field_pool_acquire_field(pool)
    copy(field_clone, field)
    return field_clone
}

face_view_clone :: proc(faces: Face_View, pool: ^Field_View_Pool) -> Face_View {
    faces_clone := field_pool_acquire_faces(pool)
    copy(faces_clone, faces)
    return faces_clone
}

view_clone :: proc {
    field_view_clone,
    face_view_clone,
}

DEFAULT_ALLOCATED_FIELDS :: 1
field_pool_init :: proc(pool: ^Field_View_Pool, mesh: Mesh, backing := context.allocator, amount := DEFAULT_ALLOCATED_FIELDS) {
    pool.backing_allocator = backing
    context.allocator = backing

    pool.num_cells = len(mesh.cells)
    pool.num_faces = len(mesh.faces)

    pool.free_fields = make([dynamic]Field_View, 0, DEFAULT_ALLOCATED_FIELDS)
    pool.free_faces = make([dynamic]Face_View, 0, DEFAULT_ALLOCATED_FIELDS)
    for _ in 0..<amount {
        append(&pool.free_fields, make(Field_View, pool.num_cells))
        append(&pool.free_faces, make(Face_View, pool.num_faces))
    }
}

field_pool_destroy :: proc(pool: ^Field_View_Pool) {
    for f in pool.free_fields{
        delete(f, pool.backing_allocator)
    }
    for f in pool.free_faces{
        delete(f, pool.backing_allocator)
    }
    delete(pool.free_fields)
    delete(pool.free_faces)
}

field_pool_acquire_field :: proc(pool: ^Field_View_Pool) -> Field_View {
   return pop_safe(&pool.free_fields) or_else make(Field_View, pool.num_cells, pool.backing_allocator)
}

field_pool_acquire_faces :: proc(pool: ^Field_View_Pool) -> Face_View {
    return pop_safe(&pool.free_faces) or_else  make(Face_View, pool.num_faces, pool.backing_allocator)
}

field_pool_acquire_vfield :: proc(pool: ^Field_View_Pool) -> Vector_Field_View {
    return {field_pool_acquire_field(pool), field_pool_acquire_field(pool)}
}

field_pool_return_fields :: proc(pool: ^Field_View_Pool, fields: ..Field_View) {
    for &field in fields {
        slice.zero(field)
        append(&pool.free_fields, field)
    }
}

field_pool_return_faces :: proc(pool: ^Field_View_Pool, faces: ..Face_View) {
    for face_field in faces {
        slice.zero(face_field)
        append(&pool.free_faces, face_field)
    }
}

field_pool_return :: proc { field_pool_return_fields, field_pool_return_faces }