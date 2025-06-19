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

import "core:math/linalg"
import "core:log"

// -- Explicit Operators --

// These will always recompute the value, in most cases it's better to use field_** functions

compute_faces_interpolate :: proc(mesh: Mesh, field: Field, faces_out: Face_View) {
    for face, face_id in mesh.faces {
        primary := mesh.cells[face.primary]
        distance_primary := (face.position - primary.position)
        if secondary_idx, has_secondary := face.secondary.?; has_secondary {
            secondary := mesh.cells[secondary_idx]
            distance_secondary := (secondary.position - face.position)
            weight := linalg.dot(face.normal, distance_secondary) / linalg.dot(face.normal, distance_secondary + distance_primary)
            faces_out[face_id] = (weight * field.data[face.primary]) + ((1 - weight) * field.data[secondary_idx])
        } else {
            bc := field.boundary_conditions[face.boundary_tag.?] or_else panic("Face has no neighbour, but does not have a boundary set.")
            switch bc.kind {
                case .Dirichlet: faces_out[face_id] = bc.value
                case .Neumann: faces_out[face_id] = field.data[face.primary] + bc.value * linalg.dot(distance_primary, face.normal)
            }
        }
    }
}

compute_flux :: proc(mesh: Mesh,  vf_faces: [2]Face_View,  flux_out: Face_View) {
    for &flux, face_id in flux_out {
        flux = linalg.dot(Vector{vf_faces.x[face_id], vf_faces.y[face_id]},  mesh.faces[face_id].normal * mesh.faces[face_id].area)
    }
}

compute_divergence :: proc(mesh: Mesh, flux: Face_View, field_out: Field_View) {
    for cell, cell_id in mesh.cells {
        flux_sum := 0.0
        for face_id in cell.faces {
            face_geometry := mesh.faces[face_id]
            is_owner := true if cell_id == face_geometry.primary else false
            flux_sum += flux[face_id] if is_owner else -flux[face_id]
        }
        field_out[cell_id] = flux_sum / cell.volume
    }
}

compute_gradient :: proc(mesh: Mesh, faces: Face_View, gradient_out: Vector_Field_View) {
    for cell, cell_id in mesh.cells {
        for face_id in cell.faces {
            face_geometry := mesh.faces[face_id]
            is_owner := true if cell_id == face_geometry.primary else false
            normal := face_geometry.normal if is_owner else -face_geometry.normal
            face_area_vector := normal * face_geometry.area
            gradient_out.x[cell_id] += faces[face_id] * face_area_vector.x / cell.volume
            gradient_out.y[cell_id] += faces[face_id] * face_area_vector.y / cell.volume
        }
    }
}

// -- Implicit Operators --

// no orthogonal correction, you can always use the corrected version even on an orthogonal mesh
add_laplacian_no_correction :: proc(lsb: ^Linear_System_Builder, mesh: Mesh, field: Field, cell_id: CellId, multiplier: f64) {
    for face_id in mesh.cells[cell_id].faces {
        ctx := get_face_context(mesh, cell_id, face_id)
        cell := mesh.cells[cell_id]
        face_area_magnitude := linalg.length(ctx.normal * ctx.face.area)
        if neighbour_id, has_neighbour := ctx.neighbour.?; has_neighbour {
            coeff := (1 / distance_cell(cell, mesh.cells[neighbour_id])) * face_area_magnitude * multiplier * f64(lsb.mode)
            lsb.row_builder[cell_id] += -coeff
            lsb.row_builder[neighbour_id] += coeff
        }else {
            bc := field.boundary_conditions[ctx.face.boundary_tag.?] or_else panic("Face has no neighbour, but does not have a boundary set.")
            switch bc.kind {
                case .Dirichlet:
                    coeff := (1 / distance_face(cell, ctx.face)) * face_area_magnitude * multiplier * f64(lsb.mode)
                    lsb.source_terms[cell_id] -= coeff * bc.value
                    lsb.row_builder[cell_id] += -coeff
                case .Neumann:
                    lsb.source_terms[cell_id] -= face_area_magnitude * bc.value * multiplier * f64(lsb.mode)
            }
        }
    }
}

//laplacian with corrections for a non-orthogonal mesh (vector from cell centroids is not colinear with face normal)
add_laplacian :: proc(lsb: ^Linear_System_Builder, mesh: Mesh, field: ^Field, cell_id: CellId, multiplier: f64) {
    grad := field_gradient(mesh, field)
    for face_id in mesh.cells[cell_id].faces {
        ctx := get_face_context(mesh, cell_id, face_id)
        cell := mesh.cells[cell_id]
        face_area_magnitude := linalg.length(ctx.normal * ctx.face.area)
        if neighbour_id, has_neighbour := ctx.neighbour.?; has_neighbour {
            ctx.face.normal = ctx.normal
            neighbour := mesh.cells[neighbour_id]
            cd := (1.0/distance_cell_ortho(cell, neighbour, ctx.face))
            coeff := cd * face_area_magnitude * multiplier * f64(lsb.mode)
            lsb.row_builder[cell_id] += -coeff
            lsb.row_builder[neighbour_id] += coeff

            //ortho correction term
            d := neighbour.position - cell.position
            d_hat := d / linalg.length(d)
            corr := ctx.normal - d_hat
            distance_primary := (ctx.face.position - cell.position)
            distance_secondary := (neighbour.position - ctx.face.position)
            weight := linalg.dot(ctx.normal, distance_secondary) / linalg.dot(ctx.normal, distance_secondary + distance_primary)
            gradxf := (weight * grad.x[cell_id]) + ((1 - weight) * grad.x[neighbour_id])
            gradyf := (weight * grad.y[cell_id]) + ((1 - weight) * grad.y[neighbour_id])
            lsb.source_terms[cell_id] -= linalg.dot(corr, [2]f64{gradxf, gradyf}) * face_area_magnitude * multiplier * f64(lsb.mode)
        }else {
            bc := field.boundary_conditions[ctx.face.boundary_tag.?] or_else panic("Face has no neighbour, but does not have a boundary set.")
            switch bc.kind {
                case .Dirichlet:
                    coeff := (1 / distance_face(cell, ctx.face)) * face_area_magnitude * multiplier * f64(lsb.mode)
                    lsb.source_terms[cell_id] -= coeff * bc.value
                    lsb.row_builder[cell_id] += -coeff
                case .Neumann:
                    lsb.source_terms[cell_id] -= face_area_magnitude * bc.value * multiplier * f64(lsb.mode)
            }
        }
    }
}


//upwind scheme advection
//TODO: more accurate scheme, like linear upwind or QUICK
add_advection :: proc(lsb: ^Linear_System_Builder, mesh: Mesh, field: Field, volume_flux: Face_View, cell_id: CellId) {
    for face_id in mesh.cells[cell_id].faces {
        ctx := get_face_context(mesh, cell_id, face_id)
        flux := volume_flux[face_id] if ctx.is_primary else -volume_flux[face_id]
        if neighbour_id, has_neighbour := ctx.neighbour.?; has_neighbour {
            lsb.row_builder[cell_id if flux >= 0 else neighbour_id] += flux * f64(lsb.mode)
        }else{
            bc := field.boundary_conditions[ctx.face.boundary_tag.?] or_else log.panic("Face has no neighbour, but does not have a boundary set.")
            switch bc.kind {
                case .Dirichlet:
                    lsb.source_terms[cell_id] -= flux * bc.value  * f64(lsb.mode)
                case .Neumann:
                    lsb.row_builder[cell_id] += flux * f64(lsb.mode)
                    lsb.source_terms[cell_id] -= bc.value * flux * distance_face(mesh.cells[cell_id], ctx.face) * f64(lsb.mode)
            }
        }
    }
}

//1st order fully implicit euler, stable but lacking in accuracy
//TODO: higher order schemes
add_time :: proc(lsb: ^Linear_System_Builder, mesh: Mesh, field: Field, cell_id: CellId, dt: f64) {
    lsb.row_builder[cell_id] += f64(lsb.mode) * mesh.cells[cell_id].volume / dt
    lsb.source_terms[cell_id] += f64(lsb.mode) * field.data[cell_id] * (mesh.cells[cell_id].volume / dt)
}

// -- Utilities --

@(private="file")
FaceContext :: struct {
    face: Face,
    is_primary: bool,
    normal: Vector,
    neighbour: Maybe(CellId),
}

@(private="file")
get_face_context :: proc(mesh: Mesh, cell_id: CellId, face_id: FaceId) -> FaceContext {
    face := mesh.faces[face_id]
    is_primary := cell_id == face.primary
    normal := face.normal if is_primary else -face.normal
    neighbour: Maybe(CellId) = face.secondary if is_primary else face.primary
    return FaceContext{face, is_primary, normal, neighbour}
}