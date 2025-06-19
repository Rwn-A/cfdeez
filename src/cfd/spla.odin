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

import "core:slice"
import "core:log"
import "base:runtime"

// Warning: This code assumes the linear system is built in order, that is coefficents for cell 0 are added, then cell 1 and so on.

Linear_System_Builder :: struct {
    row_builder: map[CellId]f64,
    current_cell: CellId,
    mode: enum {Add = 1, Sub = -1},

    values: [dynamic]f64,
    columns: [dynamic]int,
    row_indices: []int,

    allocator: runtime.Allocator,

    source_terms: []f64,
}

NUM_NEIGHBOURS_GUESS :: 5
//Note: A natural choice for this allocator would be the main case arena
linear_system_builder_new :: proc(mesh: Mesh, allocator := context.allocator) -> Linear_System_Builder {
    context.allocator = allocator
    return {
        row_builder = make(map[CellId]f64),
        current_cell = 0,
        mode = .Add,
        values = make([dynamic]f64, 0, len(mesh.cells) * NUM_NEIGHBOURS_GUESS),
        columns = make([dynamic]int, 0, len(mesh.cells) * NUM_NEIGHBOURS_GUESS),
        row_indices = make([]int, len(mesh.cells) + 1),
        source_terms = make([]f64, len(mesh.cells)),
        allocator = allocator,
    }
}

linear_system_builder_destroy :: proc(lsb: ^Linear_System_Builder) {
    delete(lsb.row_builder)
    delete(lsb.values)
    delete(lsb.columns)
    delete(lsb.row_indices, lsb.allocator)
    delete(lsb.source_terms, lsb.allocator)
}

//Warning: This does NOT clone any memory, if builder is reset, system is no longer valid.
linear_system_builder_assemble :: proc(lsb: ^Linear_System_Builder) -> Linear_System {
    return {
        M = {
            values = lsb.values[:],
            row_indices = lsb.row_indices,
            columns = lsb.columns[:],
        },
        b = lsb.source_terms,
   }
}

linear_system_builder_reset :: proc(lsb: ^Linear_System_Builder) {
    clear(&lsb.values)
    clear(&lsb.columns)
    lsb.current_cell = 0
    linear_system_builder_reset_cell(lsb)
    slice.zero(lsb.row_indices)
    slice.zero(lsb.source_terms)
}

linear_system_builder_start_cell :: proc(lsb: ^Linear_System_Builder, cell_id: CellId) {
    assert(len(lsb.row_builder) == 0, "Likely unclosed cell equation in linear system builder.")
    lsb.current_cell = cell_id
    lsb.mode = .Add
}

linear_system_builder_end_cell :: proc(lsb: ^Linear_System_Builder) {
    lsb.row_indices[lsb.current_cell] = len(lsb.values)
    defer lsb.row_indices[lsb.current_cell + 1] = len(lsb.values)
    for col, val in lsb.row_builder {
        append(&lsb.values, val)
        append(&lsb.columns, col)
    }
    linear_system_builder_reset_cell(lsb)
}

//Note: The source is multiplied by the cell volume on insertion, this makes it easier to add in fields without first having to scale them.
linear_system_builder_add_source :: proc(lsb: ^Linear_System_Builder, mesh: Mesh, source: []f64) {
    for _,  i in source {
        lsb.source_terms[i] += source[i] * f64(lsb.mode) * mesh.cells[i].volume
    }
}

@(private="file")
linear_system_builder_reset_cell :: proc(lsb: ^Linear_System_Builder) {
    clear(&lsb.row_builder)
    lsb.mode = .Add
}

Linear_System :: struct {
    M: Matrix,
    b: []f64,
}

linear_system_get_diagonal :: proc(ls: Linear_System, mesh: Mesh, out: []f64)  {
    for _, cell_id in mesh.cells {
        out[cell_id] = matrix_get_coeff(ls.M, cell_id, cell_id)
    }
}

//this is destructive, the usage of this is mostly for calculating the H vector used in PISO and SIMPLE
//since this is computed after the momentum predictor is solved, the destructive nature is fine.
linear_system_clear_diagonal :: proc(ls: ^Linear_System, mesh: Mesh) {
    for _, cell_id in mesh.cells {
        for i in ls.M.row_indices[cell_id]..<ls.M.row_indices[cell_id + 1] {
            if ls.M.columns[i] == cell_id do ls.M.values[i] = 0
        }
    }
}

//csr format
Matrix :: struct {
    row_indices: []int,
    columns: []int,
    values: []f64,
}

matrix_multiply_vector :: proc(M: Matrix, x: []f64, out: []f64) {
    for row in 0..<len(M.row_indices)-1 {
        for i in M.row_indices[row]..<M.row_indices[row + 1] {
            col := M.columns[i]
            out[row] += M.values[i] * x[col]
        }
    }
}

matrix_get_coeff :: proc(M: Matrix, row: int, column: int) -> f64 {
    for i in M.row_indices[row]..<M.row_indices[row + 1] {
        if M.columns[i] == column do return M.values[i]
    }
    return 0
}

Linsolve_Error :: enum {
    None,
    Zero_Diagonal,
    Did_Not_Converge,
}

linear_system_solve_gs :: proc(ls: Linear_System, out_field: ^Field, max_iterations := 1000, tolerance := 1e-5) -> Linsolve_Error {
    out := cast([]f64)out_field.data
    field_dirty(out_field)

    assert(len(ls.b) == len(out))
    assert(len(ls.M.row_indices) == len(ls.b) + 1)
    n := len(ls.b)

    for iter := 0; iter < max_iterations; iter += 1 {
        max_diff := 0.0
        for i in 0..<n {

            a_ii := matrix_get_coeff(ls.M, i, i)

            if abs(a_ii) < 1e-10 do return .Zero_Diagonal

            sum := 0.0
            for j in ls.M.row_indices[i]..<ls.M.row_indices[i+1] {
                col := ls.M.columns[j]
                val := ls.M.values[j]

                if col != i {
                    sum += val * out[col]
                }
            }

            old_value := out[i]
            out[i] = (ls.b[i] - sum) / a_ii

            diff := abs(out[i] - old_value)
            if diff > max_diff {
                max_diff = diff
            }
        }

        if iter % (max_iterations / 10) == 0 do log.debugf("Current residual: %e", max_diff)

        if max_diff < tolerance {
            log.debugf("Linear solve converged after %d iterations", iter)
            return .None
        }
    }

    return .Did_Not_Converge
}