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
Field_Data_Pool :: struct {
	backing:              mem.Allocator,
	num_faces, num_cells: int,
	free_cell_data:       [dynamic]Cell_Data,
	free_face_data:       [dynamic]Face_Data,
}

field_pool_init :: proc(pool: ^Field_Data_Pool, mesh: Mesh, backing := context.allocator) {
	DEFAULT_ALLOCATED_FIELDS :: 8

	pool.backing = backing
	context.allocator = backing

	pool.num_cells = len(mesh.cells)
	pool.num_faces = len(mesh.faces)

	pool.free_cell_data = make([dynamic]Cell_Data, 0, DEFAULT_ALLOCATED_FIELDS)
	pool.free_face_data = make([dynamic]Face_Data, 0, DEFAULT_ALLOCATED_FIELDS)

	for _ in 0 ..< DEFAULT_ALLOCATED_FIELDS {
		append(&pool.free_cell_data, make(Cell_Data, pool.num_cells))
		append(&pool.free_face_data, make(Face_Data, pool.num_faces))
	}
}


field_pool_destroy :: proc(pool: ^Field_Data_Pool) {
	for f in pool.free_cell_data {
		delete(f, pool.backing)
	}
	for f in pool.free_face_data {
		delete(f, pool.backing)
	}
	delete(pool.free_cell_data)
	delete(pool.free_face_data)
}

field_pool_cell_data :: proc(pool: ^Field_Data_Pool) -> Cell_Data {
	return pop_safe(&pool.free_cell_data) or_else make(Cell_Data, pool.num_cells, pool.backing)
}

field_pool_face_data :: proc(pool: ^Field_Data_Pool) -> Face_Data {
	return pop_safe(&pool.free_face_data) or_else make(Face_Data, pool.num_faces, pool.backing)
}

field_pool_vector_cell_data :: proc(pool: ^Field_Data_Pool) -> [2]Cell_Data {
	return {field_pool_cell_data(pool), field_pool_cell_data(pool)}
}

field_pool_return_cell_data :: proc(pool: ^Field_Data_Pool, items: ..Cell_Data) {
	for cell_data in items {
		slice.zero(cell_data)
		append(&pool.free_cell_data, cell_data)
	}
}

field_pool_return_face_data :: proc(pool: ^Field_Data_Pool, items: ..Face_Data) {
	for face_data in items {
		slice.zero(face_data)
		append(&pool.free_face_data, face_data)
	}
}

field_pool_return :: proc {
	field_pool_return_cell_data,
	field_pool_return_face_data,
}
