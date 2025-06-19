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

package cfd_io

import "../../cfd"

import "core:bufio"
import "core:io"
import "core:log"
import "core:mem/virtual"
import "core:math/linalg"
import "core:strings"
import "core:strconv"

//TODO: confirm performance for large meshes.

@(private="file")
Msh_Element :: struct {
    id: int,
    type: Element_Type,
    tags: []int,
    nodes: []int,
}

@(private="file")
Msh_Node :: struct {
    id: int,
    position: linalg.Vector3f64,
}

@(private="file")
Msh_Group :: struct {
    dimension: Dimension,
    id: int,
    name: string,
}

@(private="file")
Dimension :: enum { Line, Surface, Volume }

@(private="file")
Element_Type :: enum {
    MSH_LIN_2    = 1,
    MSH_TRI_3    = 2,
    MSH_QUA_4    = 3,
    MSH_TET_4    = 4,
    MSH_HEX_8    = 5,
    MSH_PRI_6    = 6,
    MSH_PYR_5    = 7,
    MSH_LIN_3    = 8,
    MSH_TRI_6    = 9,
    MSH_QUA_9    = 10,
    MSH_TET_10   = 11,
    UNSUPPORTED,
}

@(private="file")
parse_element :: proc(line: string, ln: int) -> (Msh_Element, bool) {
    context.allocator = context.temp_allocator
    components := strings.split(line, " ")
    if len(components) < 5 {
        log.errorf("Failed to read mesh file at line %d, expected atleast 5 entries in element got %s.", ln, line)
        return {}, false
    }

    id := strconv.atoi(components[0]) - 1

    type_num := strconv.atoi(components[1])
    type: Element_Type = Element_Type(type_num) if type_num < int(Element_Type.UNSUPPORTED) else .UNSUPPORTED

    num_tags := strconv.atoi(components[2])
    tags := make([]int, num_tags)
    for &t, i in tags do t = strconv.atoi(components[3 + i])

    nodes := make([]int, len(components) - (3 + num_tags))
    for &n, i in nodes do n = strconv.atoi(components[3 + num_tags + i]) - 1

    return {id, type, tags, nodes}, true
}

@(private="file")
parse_node :: proc(line: string, ln: int) -> (Msh_Node, bool) {
    context.allocator = context.temp_allocator
    components := strings.split(line, " ")
    if len(components) != 4 {
        log.errorf("Failed to read mesh file at line %d, expectected 4 entries in node got %s.", ln, line)
        return {}, false
    }
    id := strconv.atoi(components[0]) - 1
    position := linalg.Vector3f64{strconv.atof(components[1]), strconv.atof(components[2]), strconv.atof(components[3])}
    return {id, position}, true
}

@(private="file")
parse_group :: proc(line: string, ln: int) -> (Msh_Group, bool) {
    context.allocator = context.temp_allocator
    components := strings.split(line, " ")
    if len(components) != 3 {
        log.errorf("Failed to read mesh file at line %d, expected 3 entries in name got %s.", ln, line)
        return {}, false
    }
    dimension := Dimension(strconv.atoi(components[0]))
    id := strconv.atoi(components[1])
    name := strings.trim(components[2], "\"")
    name = strings.clone(name)
    return {dimension, id, name}, true
}

line_number: int

//"Clean Code" advocates proceed with caution.
from_gmsh :: proc(s: io.Stream, mesh_arena: ^virtual.Arena) -> (mesh: cfd.Mesh, ok: bool) {
    b: bufio.Reader
    bufio.reader_init(&b, s)
    defer bufio.reader_destroy(&b)

    line_number = 1

    consume_line :: proc(b: ^bufio.Reader) -> (string, bool) {
        line_number += 1
        line_slice, err := bufio.reader_read_slice(b, '\n')
        if err != nil {
            log.error("Failed to finish reading mesh file, file may have been changed while reading.")
            return "", false
        }
        line := string(line_slice)
        line = strings.trim_right_space(line)
        line = strings.trim_left_space(line)
        return line, true
    }

    scratch_arena := virtual.Arena{}
    if err := virtual.arena_init_growing(&scratch_arena); err != nil do log.panic("Failed to allocate memory!")
    defer virtual.arena_destroy(&scratch_arena)

    context.temp_allocator = virtual.arena_allocator(&scratch_arena)
    context.allocator = virtual.arena_allocator(mesh_arena)

    cell_buffer := make(#soa[dynamic]cfd.Cell)
    face_buffer := make(#soa[dynamic]cfd.Face)
    vertex_buffer := make(#soa[dynamic]cfd.Vector, context.allocator)

    boundary_names := make(map[string]cfd.Boundary_Id)

    face_set := make(map[[2]int]cfd.FaceId, context.temp_allocator)

    consume_line(&b) or_return
    version := consume_line(&b) or_return
    consume_line(&b) or_return

    if !strings.contains(version, "2.2 0 8") {
        log.errorf("Unsupported .msh format, please use version 2.2 of the format. Got %s", version)
        return {}, false
    }


    //groups
    _ = consume_line(&b) or_return
    num_groups := strconv.atoi(consume_line(&b) or_return)
    for _ in 0..<num_groups {
        group := parse_group(consume_line(&b) or_return, line_number) or_return
        boundary_names[strings.clone(group.name)] = group.id
    }
    _ = consume_line(&b) or_return

    //nodes
    _ = consume_line(&b) or_return
    num_nodes := strconv.atoi(consume_line(&b) or_return)
    for _ in 0..<num_nodes {
        node := parse_node(consume_line(&b) or_return, line_number) or_return
        append(&vertex_buffer, node.position.xy)
    }
    _ = consume_line(&b) or_return

    //elements
    _ = consume_line(&b) or_return
    non_cell_elements := 0
    num_els := strconv.atoi(consume_line(&b) or_return)
    for _ in 0..<num_els {
        element := parse_element(consume_line(&b) or_return, line_number) or_return
        #partial switch element.type {
            case .MSH_LIN_2: //(NOTE:) there is an assumption here that these elements will be encountered first, before any triangles or quads.
                v1 := element.nodes[0]; v2 := element.nodes[1]
                facekey := face_key(v1, v2)
                face_set[facekey] = len(face_buffer)
                face := face_from_nodes(vertex_buffer[v1], vertex_buffer[v2])
                face.boundary_tag = element.tags[0]
                face.primary = -1
                append(&face_buffer, face)
                non_cell_elements += 1
            case .MSH_TRI_3, .MSH_QUA_4:
                cell_face_builder := make([dynamic]cfd.FaceId)
                cell_vertices_builder := make([dynamic]cfd.VertexId)
                for node_idx, i in element.nodes {
                    node := vertex_buffer[node_idx]
                    next_node_idx := element.nodes[(i + 1) % len(element.nodes)]
                    next_node := vertex_buffer[next_node_idx]

                    facekey := face_key(node_idx, next_node_idx)
                    defer {
                        append(&cell_face_builder, face_set[facekey])
                        append(&cell_vertices_builder, node_idx)
                    }
                    if facekey in face_set {
                        if face_buffer[face_set[facekey]].primary == -1 {
                            face_buffer[face_set[facekey]].primary = element.id - non_cell_elements
                        }else{
                            face_buffer[face_set[facekey]].secondary = element.id - non_cell_elements
                        }
                        continue
                    }

                    new_face := face_from_nodes(node.xy, next_node.xy)
                    new_face.primary = element.id - non_cell_elements
                    face_set[facekey] = len(face_buffer)
                    append(&face_buffer, new_face)
                }
                new_cell: cfd.Cell
                new_cell.vertices = cell_vertices_builder[:] //we still build the buffer to use for volume and centroid
                new_cell.faces = cell_face_builder[:]
                new_cell.volume = cell_volume(vertex_buffer, cell_vertices_builder[:])
                new_cell.position = cell_center(vertex_buffer, cell_vertices_builder[:])
                append(&cell_buffer, new_cell)
            case:
                log.errorf("Failed to load mesh unsupported element type %s, the mesh must be 1st order triangles or quads", element.type)
                return {}, false
        }
    }

    return {cell_buffer[:], face_buffer[:], vertex_buffer[:], boundary_names}, true

    face_from_nodes :: proc(nodea, nodeb: cfd.Vector) -> (f: cfd.Face) {
        tangent_vector := nodeb - nodea
        f.position = (nodea + nodeb) / 2
        f.normal = linalg.normalize(cfd.Vector{tangent_vector.y, -tangent_vector.x})
        f.area = linalg.distance(nodea, nodeb)
        return f
    }

    face_key :: proc(idxa, idxb: int) -> [2]int {
        return { min(idxa,idxb), max(idxa,idxb) }
    }

    cell_volume :: proc(vertices: #soa[dynamic]cfd.Vector, vertex_ids: []cfd.VertexId) -> f64 {
        numPoints := len(vertex_ids)
        area := 0.0   // Accumulates area
        j := numPoints-1
        for i in 0..< numPoints {
            area +=  (vertices[vertex_ids[j]].x+vertices[vertex_ids[i]].x) * (vertices[vertex_ids[j]].y-vertices[vertex_ids[i]].y)
            j = i  //j is previous vertex to i
        }
        return area/-2.0
    }

    cell_center :: proc(vertices: #soa[dynamic]cfd.Vector, vertex_ids: []cfd.VertexId) -> cfd.Vector {
        new_vector: cfd.Vector
        for i in vertex_ids do new_vector += vertices[i]
        return new_vector/f64(len(vertex_ids))
    }
}