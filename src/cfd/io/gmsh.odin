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
import "core:math/linalg"
import "core:mem/virtual"
import "core:strconv"
import "core:strings"

@(private = "file")
Msh_Element :: struct {
	id:    int,
	type:  Element_Type,
	tags:  []int,
	nodes: []int,
}

@(private = "file")
Msh_Node :: struct {
	id:       int,
	position: linalg.Vector3f64,
}

@(private = "file")
Msh_Group :: struct {
	dimension: Dimension,
	id:        int,
	name:      string,
}

@(private = "file")
Dimension :: enum {
	Line,
	Surface,
	Volume,
}

@(private = "file")
Element_Type :: enum {
	MSH_LIN_2 = 1,
	MSH_TRI_3 = 2,
	MSH_QUA_4 = 3,
	MSH_TET_4 = 4,
	MSH_HEX_8 = 5,
	MSH_PRI_6 = 6,
	MSH_PYR_5 = 7,
	MSH_LIN_3 = 8,
	MSH_TRI_6 = 9,
	MSH_QUA_9 = 10,
	MSH_TET_10 = 11,
	UNSUPPORTED,
}


EXPECTED_MSH_VERSION :: "2.2 0 8"
MIN_ELEMENT_COMPONENTS :: 5
NODE_COMPONENTS :: 4
GROUP_COMPONENTS :: 3

@(private = "file")
parse_element :: proc(line: string, ln: int) -> (element: Msh_Element, ok: bool) {
	context.allocator = context.temp_allocator
	components := strings.split(line, " ")

	if len(components) < MIN_ELEMENT_COMPONENTS {
		log.errorf(
			"Failed to read mesh file at line %d, expected at least %d entries in element, got %s",
			ln,
			MIN_ELEMENT_COMPONENTS,
			line,
		)
		return {}, false
	}

	element.id = strconv.atoi(components[0]) - 1

	type_num := strconv.atoi(components[1])
	element.type = Element_Type(type_num) if type_num < int(Element_Type.UNSUPPORTED) else .UNSUPPORTED

	num_tags := strconv.atoi(components[2])
	element.tags = make([]int, num_tags)
	for &tag, i in element.tags {
		tag = strconv.atoi(components[3 + i])
	}

	num_nodes := len(components) - (3 + num_tags)
	element.nodes = make([]int, num_nodes)
	for &node, i in element.nodes {
		node = strconv.atoi(components[3 + num_tags + i]) - 1
	}

	return element, true
}

@(private = "file")
parse_node :: proc(line: string, ln: int) -> (node: Msh_Node, ok: bool) {
	context.allocator = context.temp_allocator
	components := strings.split(line, " ")

	if len(components) != NODE_COMPONENTS {
		log.errorf("Failed to read mesh file at line %d, expected %d entries in node, got %s", ln, NODE_COMPONENTS, line)
		return {}, false
	}

	node.id = strconv.atoi(components[0]) - 1
	node.position = {strconv.atof(components[1]), strconv.atof(components[2]), strconv.atof(components[3])}

	return node, true
}

@(private = "file")
parse_group :: proc(line: string, ln: int) -> (group: Msh_Group, ok: bool) {
	context.allocator = context.temp_allocator
	components := strings.split(line, " ")

	if len(components) != GROUP_COMPONENTS {
		log.errorf("Failed to read mesh file at line %d, expected %d entries in group, got %s", ln, GROUP_COMPONENTS, line)
		return {}, false
	}

	group.dimension = Dimension(strconv.atoi(components[0]))
	group.id = strconv.atoi(components[1])
	group.name = strings.clone(strings.trim(components[2], "\""))

	return group, true
}

@(private = "file")
Reader_Context :: struct {
	reader:      ^bufio.Reader,
	line_number: int,
}

@(private = "file")
consume_line :: proc(ctx: ^Reader_Context) -> (line: string, ok: bool) {
	ctx.line_number += 1
	line_slice, err := bufio.reader_read_slice(ctx.reader, '\n')
	if err != nil {
		log.error("Failed to finish reading mesh file, file may have been changed while reading.")
		return "", false
	}

	line = string(line_slice)
	line = strings.trim_space(line)
	return line, true
}

@(private = "file")
validate_version :: proc(version: string) -> bool {
	return strings.contains(version, EXPECTED_MSH_VERSION)
}

@(private = "file")
face_key :: proc(idx_a, idx_b: int) -> [2]int {
	return {min(idx_a, idx_b), max(idx_a, idx_b)}
}

@(private = "file")
process_line_element :: proc(
	element: Msh_Element,
	face_buffer: ^#soa[dynamic]cfd.Face,
	face_set: ^map[[2]int]cfd.FaceId,
	vertex_buffer: #soa[]cfd.Vector,
	non_cell_elements: ^int,
) {
	v1, v2 := element.nodes[0], element.nodes[1]
	facekey := face_key(v1, v2)
	face_set[facekey] = len(face_buffer)

	face: cfd.Face
	face.position, face.normal, face.area = cfd.face_geometry_from_vertices(vertex_buffer, {v1, v2})
	face.secondary = cfd.BoundaryId(element.tags[0])
	face.primary = -1

	append(face_buffer, face)
	non_cell_elements^ += 1
}

@(private = "file")
process_cell_element :: proc(
	element: Msh_Element,
	cell_buffer: ^#soa[dynamic]cfd.Cell,
	face_buffer: ^#soa[dynamic]cfd.Face,
	face_set: ^map[[2]int]cfd.FaceId,
	vertex_buffer: #soa[]cfd.Vector,
	non_cell_elements: int,
) {
	cell_face_builder := make([dynamic]cfd.FaceId)
	cell_vertices_builder := make([dynamic]cfd.VertexId)

	for node_idx, i in element.nodes {
		next_node_idx := element.nodes[(i + 1) % len(element.nodes)]
		facekey := face_key(node_idx, next_node_idx)

		defer {
			append(&cell_face_builder, face_set[facekey])
			append(&cell_vertices_builder, node_idx)
		}

		if facekey not_in face_set {
			face: cfd.Face
			face.position, face.normal, face.area = cfd.face_geometry_from_vertices(
				vertex_buffer,
				{node_idx, next_node_idx},
			)
			face.primary = element.id - non_cell_elements
			face_set[facekey] = len(face_buffer)
			append(face_buffer, face)
			continue
		}

		existing_face := &face_buffer[face_set[facekey]]
		if existing_face.primary == -1 {
			existing_face.primary = element.id - non_cell_elements
		} else {
			existing_face.secondary = cfd.CellId(element.id - non_cell_elements)
		}
	}

	new_cell: cfd.Cell
	new_cell.vertices = cell_vertices_builder[:]
	new_cell.faces = cell_face_builder[:]
	new_cell.volume = cfd.cell_volume_from_vertices(vertex_buffer, cell_vertices_builder[:])
	new_cell.position = cfd.cell_center_from_vertices(vertex_buffer, cell_vertices_builder[:])

	append(cell_buffer, new_cell)
}

from_gmsh :: proc(
	s: io.Stream,
	mesh_arena: ^virtual.Arena,
) -> (
	mesh: cfd.Mesh,
	names: map[string]cfd.BoundaryId,
	ok: bool,
) {

	b: bufio.Reader
	bufio.reader_init(&b, s)
	defer bufio.reader_destroy(&b)

	ctx := Reader_Context {
		reader      = &b,
		line_number = 1,
	}


	scratch_arena := virtual.Arena{}
	if err := virtual.arena_init_growing(&scratch_arena); err != nil {
		log.panic("Failed to allocate memory!")
	}
	defer virtual.arena_destroy(&scratch_arena)

	context.temp_allocator = virtual.arena_allocator(&scratch_arena)
	context.allocator = virtual.arena_allocator(mesh_arena)

	cell_buffer := make(#soa[dynamic]cfd.Cell)
	face_buffer := make(#soa[dynamic]cfd.Face)
	vertex_buffer := make(#soa[dynamic]cfd.Vector)
	boundary_names := make(map[string]cfd.BoundaryId)
	face_set := make(map[[2]int]cfd.FaceId, context.temp_allocator)

	consume_line(&ctx) or_return
	version := consume_line(&ctx) or_return
	consume_line(&ctx) or_return

	if !validate_version(version) {
		log.errorf("Unsupported .msh format, please use version %s of the format. Got %s", EXPECTED_MSH_VERSION, version)
		return {}, {}, false
	}

	consume_line(&ctx) or_return // $PhysicalNames
	num_groups := strconv.atoi(consume_line(&ctx) or_return)

	for _ in 0 ..< num_groups {
		group := parse_group(consume_line(&ctx) or_return, ctx.line_number) or_return
		boundary_names[strings.clone(group.name)] = cfd.BoundaryId(group.id)
	}
	consume_line(&ctx) or_return // $EndPhysicalNames

	// Parse nodes
	consume_line(&ctx) or_return // $Nodes
	num_nodes := strconv.atoi(consume_line(&ctx) or_return)

	for _ in 0 ..< num_nodes {
		node := parse_node(consume_line(&ctx) or_return, ctx.line_number) or_return
		append(&vertex_buffer, node.position.xy)
	}
	consume_line(&ctx) or_return // $EndNodes

	consume_line(&ctx) or_return // $Elements
	num_elements := strconv.atoi(consume_line(&ctx) or_return)
	non_cell_elements := 0

	for _ in 0 ..< num_elements {
		element := parse_element(consume_line(&ctx) or_return, ctx.line_number) or_return

		#partial switch element.type {
		case .MSH_LIN_2:
			process_line_element(element, &face_buffer, &face_set, vertex_buffer[:], &non_cell_elements)
		case .MSH_TRI_3, .MSH_QUA_4:
			process_cell_element(element, &cell_buffer, &face_buffer, &face_set, vertex_buffer[:], non_cell_elements)
		case:
			log.errorf(
				"Failed to load mesh: unsupported element type %v. The mesh must contain 1st order triangles or quads",
				element.type,
			)
			return {}, {}, false
		}
	}

	return {cell_buffer[:], face_buffer[:], vertex_buffer[:]}, boundary_names, true
}
