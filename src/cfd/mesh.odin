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

Boundary_Id :: int

Mesh :: struct {
    cells: #soa[]Cell,
    faces: #soa[]Face,
    vertices: #soa[]Vector,
    boundary_names: map[string]Boundary_Id,
}

FaceId :: int
CellId :: int
VertexId :: int

Cell :: struct {
    position: Vector,
    faces: []FaceId,
    volume: f64,
    vertices: []VertexId,
}

Face :: struct {
    position, normal: Vector,
    primary: CellId,
    secondary: Maybe(CellId),
    boundary_tag: Maybe(Boundary_Id),
    area: f64,
}


distance_cell :: proc(P, N: Cell) -> f64 {
    return linalg.distance(P.position, N.position)
}

//distance between cells projected onto face normal. Used in orthogonal corrections
//equivalent to distance_cell in a structured mesh.
distance_cell_ortho :: proc(P, N: Cell, f: Face) -> f64 {
    return linalg.dot(N.position - P.position, f.normal)
}


distance_face :: proc(P: Cell, f: Face) -> f64 {
    d := f.normal * (linalg.dot(f.normal, f.position - P.position))
    return linalg.length(d)
}


//keeping this around in the case of changes to mesh generation

// verify_face_normals :: proc(mesh: Mesh) {
//  for cell, id in mesh.cells{
//     for f in cell.faces{
//             face := mesh.faces[f]
//             if face.primary != id do continue
//             ctf := face.position - cell.position
//             assert(linalg.dot(ctf, face.normal) >= 0)
//             assert(linalg.length(face.normal) == 1)
//         }
//     }
// }