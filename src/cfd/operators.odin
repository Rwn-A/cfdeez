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

import "core:log"
import "core:math/linalg"
import "core:slice"

// Math References: https://doc.cfd.direct/notes/cfd-general-principles
// Implementation entirely our own.

// -- Explicit operations, requires only known field values --

compute_value_across_face :: proc(mesh: Mesh, val1, val2: f64, cid1, cid2: CellId, face: Face) -> f64 {
	c1 := mesh.cells[cid1]; c2 := mesh.cells[cid2]
	face := face_ensure_outward(c1, face)
	d1 := cell_face_delta(c1, face)
	d2 := -cell_face_delta(c2, face)
	weight := linalg.dot(face.normal, d2) / linalg.dot(face.normal, d2 + d1)
	return (weight * val1) + ((1 - weight) * val2)
}

compute_faces_interpolate :: proc(mesh: Mesh, field: Scalar_Field, faces_out: Face_Data) {
	for face, face_id in mesh.faces {
		distance_primary := cell_face_delta(mesh.cells[face.primary], face)

		switch id in face.secondary {
		case BoundaryId:
			bc := field.bnds[id] or_else log.panicf("Boundary id %d was not set on field.", id)
			switch bc.kind {
			case .Dirichlet:
				faces_out[face_id] = bc.value
			case .Neumann:
				faces_out[face_id] = field.data[face.primary] + bc.value * linalg.dot(distance_primary, face.normal)
			}
		case CellId:
			val := compute_value_across_face(mesh, field.data[face.primary], field.data[id], face.primary, id, face)
			faces_out[face_id] = val
		}
	}
}

compute_flux :: proc(mesh: Mesh, vf_faces: [VF_COMPONENTS]Face_Data, flux_out: Face_Data) {
	for &flux, face_id in flux_out {
		flux = linalg.dot(Vector{vf_faces.x[face_id], vf_faces.y[face_id]}, face_area_v(mesh.faces[face_id]))
	}
}

compute_divergence :: proc(mesh: Mesh, flux: Face_Data, field_out: []f64) {
	for cell, cell_id in mesh.cells {
		flux_sum := 0.0
		for face_id in cell.faces {
			flux_sum += flux[face_id] if cell_is_owner(cell_id, mesh.faces[face_id]) else -flux[face_id]
		}
		field_out[cell_id] = flux_sum / cell.volume
	}
}

compute_gradient :: proc(mesh: Mesh, faces: Face_Data, field_data: []f64, gradient_out: [VF_COMPONENTS][]f64) {
	for cell, cell_id in mesh.cells {
		gradient_out.x[cell_id] = 0;gradient_out.y[cell_id] = 0
		for face_id in cell.faces {
			face_geometry := face_ensure_outward(cell, mesh.faces[face_id])
			face_area_vector := face_area_v(face_geometry)
			gradient_out.x[cell_id] += faces[face_id] * face_area_vector.x / cell.volume
			gradient_out.y[cell_id] += faces[face_id] * face_area_vector.y / cell.volume
		}
	}
	apply_gradient_limiting(mesh, field_data, gradient_out)
}


@(private = "file")
apply_gradient_limiting :: proc(mesh: Mesh, field_data: []f64, gradient_out: [VF_COMPONENTS][]f64) {
	for cell, cell_id in mesh.cells {
		cell_value := field_data[cell_id]
		min_neighbor := cell_value
		max_neighbor := cell_value

		for face_id in cell.faces {
			face := mesh.faces[face_id]
			neighbour_id := cell_neighbour(cell_id, face) or_continue
			neighbour_value := field_data[neighbour_id]

			min_neighbor = min(min_neighbor, neighbour_value)
			max_neighbor = max(max_neighbor, neighbour_value)
		}

		limiter := 1.0

		for face_id in cell.faces {
			face := mesh.faces[face_id]
			neighbour_id := cell_neighbour(cell_id, face) or_continue

			d := cell_delta(mesh.cells[cell_id], mesh.cells[neighbour_id])
			extrapolated_value := cell_value + gradient_out.x[cell_id] * d.x + gradient_out.y[cell_id] * d.y
			if extrapolated_value > max_neighbor {
				face_limiter := (max_neighbor - cell_value) / (extrapolated_value - cell_value)
				limiter = min(limiter, face_limiter)
			} else if extrapolated_value < min_neighbor {
				face_limiter := (min_neighbor - cell_value) / (extrapolated_value - cell_value)
				limiter = min(limiter, face_limiter)
			}
		}
		limiter = max(0.0, limiter)
		gradient_out.x[cell_id] *= limiter
		gradient_out.y[cell_id] *= limiter
	}
}

// -- Implicit operators, to create equations for currently unknown field values --

add_laplacian :: proc(lsb: ^Linear_System_Builder, mesh: Mesh, field: ^Scalar_Field, cell_id: CellId, multiplier: f64) {
	grad := field_gradient(mesh, field)
	cell := mesh.cells[cell_id]
	for face_id in mesh.cells[cell_id].faces {
		face := face_ensure_outward(cell, mesh.faces[face_id])
		face_area_mag := linalg.length(face_area_v(face))

		if neighbour_id, has := cell_neighbour(cell_id, face); has {
			neighbour := mesh.cells[neighbour_id]
			cd := (1.0 / cell_ortho_distance(cell, neighbour, face))
			coeff := cd * face_area_mag * multiplier * f64(lsb.mode)
			lsb.row_builder[cell_id] += -coeff
			lsb.row_builder[neighbour_id] += coeff

			//ortho correction term
			d_unit := linalg.normalize(cell_delta(cell, neighbour))
			corr := face.normal - d_unit
			gradxf := compute_value_across_face(mesh, grad.x[cell_id], grad.x[neighbour_id], cell_id, neighbour_id, face)
			gradyf := compute_value_across_face(mesh, grad.y[cell_id], grad.y[neighbour_id], cell_id, neighbour_id, face)
			lsb.source_terms[cell_id] -=
				linalg.dot(corr, [2]f64{gradxf, gradyf}) * face_area_mag * multiplier * f64(lsb.mode)
		} else {
			id := face.secondary.(BoundaryId)
			bc := field.bnds[id] or_else log.panicf("Boundary id %d was not set on field.", id)
			switch bc.kind {
			case .Dirichlet:
				coeff := (1.0 / cell_face_ortho_distance(cell, face)) * face_area_mag * multiplier * f64(lsb.mode)
				lsb.source_terms[cell_id] -= coeff * bc.value
				lsb.row_builder[cell_id] += -coeff
			case .Neumann:
				lsb.source_terms[cell_id] -= face_area_mag * bc.value * multiplier * f64(lsb.mode)
			}

		}
	}
}


add_advection :: proc(
	lsb: ^Linear_System_Builder,
	mesh: Mesh,
	field: Scalar_Field,
	volume_flux: Face_Data,
	cell_id: CellId,
) {
	cell :=  mesh.cells[cell_id]
	for face_id in cell.faces{
		face := face_ensure_outward(cell, mesh.faces[face_id])
		flux := volume_flux[face_id] if cell_is_owner(cell_id, face) else -volume_flux[face_id]
		if neighbour_id, has := cell_neighbour(cell_id, face); has {
			lsb.row_builder[cell_id if flux >= 0 else neighbour_id] += flux * f64(lsb.mode)
		} else {
			id := face.secondary.(BoundaryId)
			bc := field.bnds[id] or_else log.panicf("Boundary id %d was not set on field.", id)
			switch bc.kind {
			case .Dirichlet:
				lsb.source_terms[cell_id] -= flux * bc.value * f64(lsb.mode)
			case .Neumann:
				lsb.row_builder[cell_id] += flux * f64(lsb.mode)
				lsb.source_terms[cell_id] -=
					bc.value * flux * cell_face_ortho_distance(mesh.cells[cell_id], face) * f64(lsb.mode)
			}
		}
	}
}

add_time :: proc(lsb: ^Linear_System_Builder, mesh: Mesh, field: Scalar_Field, cell_id: CellId, dt: f64) {
	lsb.row_builder[cell_id] += f64(lsb.mode) * mesh.cells[cell_id].volume / dt
	lsb.source_terms[cell_id] += f64(lsb.mode) * field.data[cell_id] * (mesh.cells[cell_id].volume / dt)
}
