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

import "core:os"
import "core:log"
import "core:fmt"
import "core:path/filepath"
import sa"core:container/small_array"
import "core:strconv"

import "../../cfd"

// TODO: IO may need to be buffered
// TODO: case name should be converted to lowercase and no spaces.

Output_Fn :: #type proc(
    mesh: cfd.Mesh,
    fields: cfd.Primary_Fields,
    passive_names: []string,
    directory, case_name: string,
    step: int
) -> bool

output_path :: proc(dir: string, case_name: string, ext: string, step: int) -> string {
    @static output_path_bufffer: [128]byte
    return fmt.bprintf(output_path_bufffer[:], "%s%c%s_%d.%s",dir,filepath.SEPARATOR, case_name, step, ext)
}

to_csv :: proc(
    mesh: cfd.Mesh,
    fields: cfd.Primary_Fields,
    passive_names: []string,
    directory, case_name: string,
    step: int
) -> bool {
    output_path := output_path(directory, case_name, "csv", step)
    file, err := os.open(output_path, os.O_CREATE | os.O_TRUNC | os.O_RDWR, 0o777)
    if err != nil {
        log.errorf("Could not create output CSV file %s", output_path)
        return false
    }
    defer os.close(file)

    fmt.fprintf(file, "x,y,u,v,p")
    for passive_name in passive_names {
        fmt.fprintf(file, ",%s", passive_name)
    }
    fmt.fprintf(file, "\n")
    for cell, i in mesh.cells {
        fmt.fprintf(file, "%f,%f,%f,%f,%f",
            cell.position.x,
            cell.position.y,
            fields.u.components.x.data[i],
            fields.u.components.y.data[i],
            fields.p.data[i],
        )
        for passive in fields.passives {
            fmt.fprintf(file, ",%f", passive.data[i])
        }
        fmt.fprintf(file, "\n")
    }
    return true
}

@(private="file")
Tag :: enum u8 {
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

@(private="file")
MAX_NESTED_TAGS :: 16

VTU_Writer :: struct {
    open: sa.Small_Array(MAX_NESTED_TAGS, Tag),
    fd: os.Handle,
}

@(private="file")
open_tag :: proc(w: ^VTU_Writer, t: Tag, attrib: [][2]string = {}) {
    sa.push(&w.open, t)
    fmt.fprintf(w.fd, "<%s", t)
    for a in attrib {
        fmt.fprintf(w.fd, " %s= \"%s\"", a[0], a[1])
    }
    fmt.fprintfln(w.fd, ">")
}

@(private="file")
close_tag :: proc(w: ^VTU_Writer) {
    fmt.fprintfln(w.fd, "</%s>", sa.pop_back(&w.open))
}

// TODO: Add Binary data instead of ascii, i think we can use vendor:zlib and core:encoding/base64
to_vtu :: proc(
    mesh: cfd.Mesh,
    fields: cfd.Primary_Fields,
    passive_names: []string,
    directory, case_name: string,
    step: int
) -> bool {
   output_path := output_path(directory, case_name, "vtu", step)
   fd, err := os.open(output_path, os.O_CREATE | os.O_TRUNC | os.O_RDWR, 0o777)
    if err != nil {
        log.errorf("Could not create output VTU file %s", output_path)
        return false
    }
    defer os.close(fd)

    w := VTU_Writer{fd = fd}

    fmt.fprintfln(fd, "<?xml version=\"1.0\"?>")

    open_tag(&w, .VTKFile, {
        {"type", "UnstructuredGrid"}, {"version", "0.1"}, {"byte_order", "LittleEndian"},
    })
    defer close_tag(&w)

    open_tag(&w, .UnstructuredGrid)
    defer close_tag(&w)

    b0: [16]u8
    b1: [16]u8
    num_points := strconv.itoa(b0[:], len(mesh.vertices))
    num_cells := strconv.itoa(b1[:], len(mesh.cells))

    open_tag(&w, .Piece, {{"NumberOfPoints", num_points}, {"NumberOfCells", num_cells}})
    defer close_tag(&w)

    {
        open_tag(&w, .Points)
        defer close_tag(&w)

        open_tag(&w, .DataArray, {
            {"type", "Float64"}, {"Name", "Points"}, {"NumberOfComponents", "3"}, {"format", "ascii"},
        })
        for v in mesh.vertices{
            fmt.fprintfln(w.fd, "%e", v.x)
            fmt.fprintfln(w.fd, "%e", v.y)
            fmt.fprintfln(w.fd, "%e", f64(0))
        }
        defer close_tag(&w)

    }

    {
        open_tag(&w, .Cells)
        defer close_tag(&w)

        open_tag(&w, .DataArray, {
            {"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"},
        })

        cell_offsets := make([dynamic]int)
        cur_offset := 0
        defer delete(cell_offsets)

        for cell in mesh.cells {
            for v in cell.vertices{
                fmt.fprintfln(w.fd, "%d", v)
                cur_offset += 1
            }
            append(&cell_offsets, cur_offset)
        }
        close_tag(&w)

        open_tag(&w, .DataArray, {
            {"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"},
        })
        for o in cell_offsets{
            fmt.fprintfln(w.fd, "%d", o)
        }
        close_tag(&w)

        open_tag(&w, .DataArray, {
            {"type", "Int64"}, {"Name", "types"}, {"format", "ascii"},
        })
        Type :: enum {Tri = 5, Quad = 9}
        for cell in mesh.cells {
            switch len(cell.vertices) {
                case 3: fmt.fprintfln(w.fd, "%d", i32(Type.Tri))
                case 4: fmt.fprintfln(w.fd, "%d", i32(Type.Quad))
                case: panic("Unsupported element in vtk output, this means the mesh import validation failed.")
            }
        }
        close_tag(&w)
    }

    {
        open_tag(&w, .CellData)
        defer close_tag(&w)

        open_tag(&w, .DataArray, {
            {"type", "Float64"}, {"Name", "velocity"}, {"NumberOfComponents", "3"}, {"format", "ascii"},
        })
        for _, i in mesh.cells{
            fmt.fprintfln(w.fd, "%e", fields.u.components.x.data[i])
            fmt.fprintfln(w.fd, "%e", fields.u.components.y.data[i])
            fmt.fprintfln(w.fd, "%e", f64(0))

        }
        close_tag(&w)

        open_tag(&w, .DataArray, {
            {"type", "Float64"}, {"Name", "pressure"}, {"format", "ascii"},
        })
        for p in fields.p.data{
            fmt.fprintfln(w.fd, "%e", p)
        }
        close_tag(&w)

        for passive, i in fields.passives {
            open_tag(&w, .DataArray, {
                {"type", "Float64"}, {"Name", passive_names[i]}, {"format", "ascii"},
            })
            defer close_tag(&w)
            for d in passive.data{
                fmt.fprintfln(w.fd, "%e", d)
            }

        }
    }

    return true
}


write_pvd_header :: proc(fd: os.Handle) -> VTU_Writer {
    w := VTU_Writer{fd = fd}
    fmt.fprintfln(fd, "<?xml version=\"1.0\"?>")
    open_tag(&w, .VTKFile, {
        {"type", "Collection"}, {"version", "0.1"}, {"byte_order", "LittleEndian"},
    })

    open_tag(&w, .Collection)

    return w
}

write_pvd_entry :: proc(w: ^VTU_Writer, filename: string, time: f64) {
    //f.write(f'    <DataSet timestep="{timestep}" group="" part="0" file="{fname}"/>\n')
    strconv_buf: [16]u8
    open_tag(w, .DataSet, {
        {"timestep", strconv.ftoa(strconv_buf[:], time, 'f', 8, 64)},
        {"group", ""},
        {"part", "0"},
        {"file", filename}
    })
    close_tag(w)
}

write_pvd_end :: proc(w: ^VTU_Writer) {
    close_tag(w) //close Collection
    close_tag(w) //close VTKFiel
}