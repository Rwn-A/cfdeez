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

Vector :: linalg.Vector2f64

Mesh :: struct {
	cells:    #soa[]Cell,
	faces:    #soa[]Face,
	vertices: #soa[]Vector, // needed strictly for output files.
}

BoundaryId :: distinct int
FaceId :: int
CellId :: int
VertexId :: int

Cell :: struct {
	position: Vector,
	faces:    []FaceId,
	vertices: []VertexId, // needed strictly for output files.
	volume:   f64,
}

Face :: struct {
	position, normal: Vector,
	primary:          CellId,
	secondary:        union #no_nil {
		CellId,
		BoundaryId,
	},
	area:             f64,
}

// -- geometry and access helpers --

cell_delta :: proc(p, n: Cell) -> Vector {
	return n.position - p.position
}

cell_face_delta :: proc(c: Cell, f: Face) -> Vector {
	return f.position - c.position
}

// distance between cells but projected onto the face normal.
// face normal is not ensured to be outward from cell p, this is signed.
cell_ortho_distance :: proc(p, n: Cell, f: Face) -> f64 {
	return linalg.dot(cell_delta(p, n), f.normal)
}

// distance between cell and face but along face normal direction, not signed.
cell_face_ortho_distance :: proc(p: Cell, f: Face) -> f64 {
	d := f.normal * (linalg.dot(f.normal, f.position - p.position))
	return linalg.length(d)
}

//returns a face with the normal pointing outward from the given cell
face_ensure_outward :: proc(cell: Cell, f: Face) -> Face {
	if linalg.dot(cell_face_delta(cell, f), f.normal) > 0 {return f}
	f := f
	f.normal *= -1
	return f
}

// returns false if cell has no neighbour
cell_neighbour :: proc(cid: CellId, f: Face) -> (CellId, bool) {
	if cid == f.primary {return f.secondary.(CellId)}
	return f.primary, true
}

cell_is_owner :: proc(cid: CellId, f: Face) -> bool {
	return cid == f.primary 
}

face_area_v :: proc(face: Face) -> Vector {
	return face.normal * face.area
}

// -- cell geometry buillers --
// Not used by simulation, handy to have for external mesh imports --
// TODO: these are correct for 2d, for 3d i'm not convinced
// implement: https://doc.cfd.direct/notes/cfd-general-principles/finite-volume-mesh

//assumes vertices in counter clockwise order.
cell_volume_from_vertices :: proc(vertices: #soa[]Vector, vertex_ids: []VertexId) -> f64 {
	volume := 0.0
	j := len(vertex_ids) - 1 // previous vertex to i
	for i in 0 ..< len(vertex_ids) {
		volume +=
			(vertices[vertex_ids[j]].x + vertices[vertex_ids[i]].x) *
			(vertices[vertex_ids[j]].y - vertices[vertex_ids[i]].y)
		j = i
	}
	return volume / -2.0
}

cell_center_from_vertices :: proc(vertices: #soa[]Vector, vertex_ids: []VertexId) -> (center: Vector) {
	for id in vertex_ids {center += vertices[id]}
	return center / cast(f64)len(vertex_ids)
}

face_geometry_from_vertices :: proc(vertices: #soa[]Vector, vertex_ids: []VertexId) -> (center, normal: Vector, area: f64) {
	assert(len(vertex_ids) == 2, "Tried to get face geometry without supplying exactly 2 vertices")
	v1 := vertices[vertex_ids[0]];v2 := vertices[vertex_ids[1]]
	tangent := v2 - v1
	center = (v1 + v2) / 2
	normal = linalg.normalize(Vector{tangent.y, -tangent.x})
	area = linalg.distance(v1, v2)
	return
}


verify_face_normals :: proc(mesh: Mesh) {
 for cell, id in mesh.cells{
    for f in cell.faces{
            face := mesh.faces[f]
            if face.primary != id do continue
            ctf := face.position - cell.position
            assert(linalg.dot(ctf, face.normal) >= 0)
            assert(linalg.length(face.normal) == 1)
        }
    }
}