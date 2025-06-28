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
import "core:fmt"
import "core:os"
import "core:bufio"
import "core:io"
import sa"core:container/small_array"
import "core:log"
import "core:strconv"
import "core:encoding/base64"
import "core:mem"

VTK_Tag :: enum u8 {
    VTKFile,
    UnstructuredGrid,
    Piece,
    Points,
    DataArray,
    Cells,
    CellData,
    Collection,
    DataSet,
}

MAX_NESTED_TAGS :: 16

VTK_Writer :: struct {
    open: sa.Small_Array(MAX_NESTED_TAGS, VTK_Tag),
    handle: os.Handle,
}

vtk_writer_init :: proc(path: string) -> (w: VTK_Writer, ok: bool) {
    fd, err := os.open(path, os.O_CREATE | os.O_TRUNC | os.O_RDWR, 0o777)
    if err != nil {
        log.errorf("Could not create file %s", path)
        return {}, false
    }
    w.handle = fd
    return w, true
}

vtk_writer_destroy :: proc(vtk_w: ^VTK_Writer) {
    os.close(vtk_w.handle)
}

vtk_write_binary_data :: proc(vtk_w: ^VTK_Writer, data: []u8) {
    data_size := u32(len(data))
    
    // have to have room for length information at start of blob
    total_size := size_of(u32) + len(data)
    binary_block := make([]u8, total_size)
    defer delete(binary_block)
    
    mem.copy_non_overlapping(raw_data(binary_block[0:4]), &data_size, size_of(u32))
    mem.copy_non_overlapping(raw_data(binary_block[4:]), raw_data(data), len(data))
    
    encoded := base64.encode(binary_block)
    defer delete(encoded)
    
    fmt.fprint(vtk_w.handle, encoded)
}

@(private="file")
slice_to_bytes :: proc(data: []$T) -> []u8 {
    return mem.slice_data_cast([]u8, data)
}

vtk_open_tag :: proc(vtk_w: ^VTK_Writer, t: VTK_Tag, attrib: [][2]string = {}) {
    sa.push(&vtk_w.open, t)
    fmt.fprintf(vtk_w.handle, "<%s", t)
    for a in attrib {
        fmt.fprintf(vtk_w.handle, " %s= \"%s\"", a[0], a[1])
    }
    fmt.fprintln(vtk_w.handle, ">")
}

vtk_close_tag :: proc(vtk_w: ^VTK_Writer) {
    fmt.fprintfln(vtk_w.handle, "</%s>", sa.pop_back(&vtk_w.open))
}

vtk_write_vtu :: proc(
    mesh: cfd.Mesh,
    vfs: []cfd.Vector_Field,
    sfs: []cfd.Scalar_Field,
    vf_names, sf_names: []string,
    directory, sim_name: string,
    step: int,
) -> (out: string, ok: bool) {
   output_path := output_path(directory, sim_name, "vtu", step)

    vtk_w := vtk_writer_init(output_path) or_return
    defer vtk_writer_destroy(&vtk_w)

    fmt.fprintln(vtk_w.handle, "<?xml version=\"1.0\"?>")

    vtk_open_tag(&vtk_w, .VTKFile, {
        {"type", "UnstructuredGrid"}, 
        {"version", "0.1"}, 
        {"byte_order", "LittleEndian"},
    })
    defer vtk_close_tag(&vtk_w)

    vtk_open_tag(&vtk_w, .UnstructuredGrid)
    defer vtk_close_tag(&vtk_w)

    b0: [24]u8
    b1: [24]u8
    num_points := strconv.itoa(b0[:], len(mesh.vertices))
    num_cells := strconv.itoa(b1[:], len(mesh.cells))

    //vertices
    vtk_open_tag(&vtk_w, .Piece, {{"NumberOfPoints", num_points}, {"NumberOfCells", num_cells}})
    defer vtk_close_tag(&vtk_w)
    {
        vtk_open_tag(&vtk_w, .Points)
        defer vtk_close_tag(&vtk_w)

        vtk_open_tag(&vtk_w, .DataArray, {
            {"type", "Float64"}, 
            {"Name", "Points"}, 
            {"NumberOfComponents", "3"}, 
            {"format", "binary"},
        })
        
        points_data := make([]f64, len(mesh.vertices) * 3)
        defer delete(points_data)
        
        for v, i in mesh.vertices {
            points_data[i*3 + 0] = v.x
            points_data[i*3 + 1] = v.y
            points_data[i*3 + 2] = 0.0
        }
        
        vtk_write_binary_data(&vtk_w, slice_to_bytes(points_data))
        vtk_close_tag(&vtk_w)
    }

    // mesh connectivity
    {
        vtk_open_tag(&vtk_w, .Cells)
        defer vtk_close_tag(&vtk_w)

        vtk_open_tag(&vtk_w, .DataArray, {
            {"type", "Int32"}, 
            {"Name", "connectivity"}, 
            {"format", "binary"},
        })

        connectivity_data := make([dynamic]i32)
        defer delete(connectivity_data)
        
        for cell in mesh.cells {
            for v in cell.vertices {
                append(&connectivity_data, i32(v))
            }
        }
        
        vtk_write_binary_data(&vtk_w, slice_to_bytes(connectivity_data[:]))
        vtk_close_tag(&vtk_w)

        vtk_open_tag(&vtk_w, .DataArray, {
            {"type", "Int32"}, 
            {"Name", "offsets"}, 
            {"format", "binary"},
        })
        
        offsets_data := make([]i32, len(mesh.cells))
        defer delete(offsets_data)
        
        cur_offset: i32 = 0
        for cell, i in mesh.cells {
            cur_offset += i32(len(cell.vertices))
            offsets_data[i] = cur_offset
        }
        
        vtk_write_binary_data(&vtk_w, slice_to_bytes(offsets_data))
        vtk_close_tag(&vtk_w)

        vtk_open_tag(&vtk_w, .DataArray, {
            {"type", "Int32"}, 
            {"Name", "types"}, 
            {"format", "binary"},
        })
        
        Type :: enum i32 {Tri = 5, Quad = 9}
        types_data := make([]i32, len(mesh.cells))
        defer delete(types_data)
        
        for cell, i in mesh.cells {
            switch len(cell.vertices) {
                case 3: types_data[i] = i32(Type.Tri)
                case 4: types_data[i] = i32(Type.Quad)
                case: panic("Unsupported element in vtk output, this means the mesh import validation failed.")
            }
        }
        
        vtk_write_binary_data(&vtk_w, slice_to_bytes(types_data))
        vtk_close_tag(&vtk_w)
    }

    // field values
    {
        vtk_open_tag(&vtk_w, .CellData)
        defer vtk_close_tag(&vtk_w)

        for vf, i in vfs {
            write_vf_binary(&vtk_w, vf, vf_names[i])
        }
        for sf, i in sfs {
            write_sf_binary(&vtk_w, sf, sf_names[i])
        }
    }

    return output_path, true

    write_vf_binary :: proc(vtk_w: ^VTK_Writer, vf: cfd.Vector_Field, name: string) {
        vtk_open_tag(vtk_w, .DataArray, {
            {"type", "Float64"}, 
            {"Name", name}, 
            {"NumberOfComponents", "3"}, 
            {"format", "binary"},
        })
        defer vtk_close_tag(vtk_w)

        vector_data := make([]f64, len(vf.components.x.data) * 3)
        defer delete(vector_data)
        
        for _, i in vf.components.x.data {
            vec := cfd.vector_field_at(vf, i)
            vector_data[i*3 + 0] = vec.x
            vector_data[i*3 + 1] = vec.y
            vector_data[i*3 + 2] = 0.0
        }
        
        vtk_write_binary_data(vtk_w, slice_to_bytes(vector_data))
    }

    write_sf_binary :: proc(vtk_w: ^VTK_Writer, sf: cfd.Scalar_Field, name: string) {
        vtk_open_tag(vtk_w, .DataArray, {
            {"type", "Float64"}, 
            {"Name", name}, 
            {"format", "binary"},
        })
        defer vtk_close_tag(vtk_w)
        
        vtk_write_binary_data(vtk_w, slice_to_bytes(sf.data))
    }
}


vtk_start_pvd :: proc(directory, sim_name: string) -> (w: VTK_Writer, ok: bool) {
    path := output_path(directory, sim_name, "pvd", 1)
    w = vtk_writer_init(path) or_return
    fmt.fprintln(w.handle, "<?xml version=\"1.0\"?>")
    vtk_open_tag(&w, .VTKFile, {
        {"type", "Collection"}, {"version", "0.1"}, {"byte_order", "LittleEndian"},
    })
    vtk_open_tag(&w, .Collection)
    return w, true
}

vtk_write_pvd_entry :: proc(vtk_w: ^VTK_Writer, filename: string, time: f64) {
    strconv_buf: [24]u8
    vtk_open_tag(vtk_w, .DataSet, {
        {"timestep", strconv.ftoa(strconv_buf[:], time, 'f', 8, 64)},
        {"group", ""},
        {"part", "0"},
        {"file", filename}
    })
    vtk_close_tag(vtk_w)
}

vtk_end_pvd :: proc(vtk_w: ^VTK_Writer) {
    vtk_close_tag(vtk_w) //close Collection
    vtk_close_tag(vtk_w) //close VTKFile
    vtk_writer_destroy(vtk_w)
}