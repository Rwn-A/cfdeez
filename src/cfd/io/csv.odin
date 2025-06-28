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
import "core:os"
import "core:log"
import "core:fmt"

write_csv :: proc(
    mesh: cfd.Mesh,
    vfs: []cfd.Vector_Field,
    sfs: []cfd.Scalar_Field,
    vf_names, sf_names: []string,
    directory, sim_name: string,
    step: int,
) -> (out: string, ok: bool) {
    output_path := output_path(directory, sim_name, "csv", step)
    file, err := os.open(output_path, os.O_CREATE | os.O_TRUNC | os.O_RDWR, 0o777)
    if err != nil {
        log.errorf("Could not create output CSV file %s", output_path)
        return "", false
    }
    defer os.close(file)


    fmt.fprintf(file, "x,y")
    for name in vf_names {
    	fmt.fprintf(file, ",%s.x", name)
    	fmt.fprintf(file, ",%s.y", name)
    }
    for name in sf_names {
    	fmt.fprintf(file, ",%", name)
    }
    fmt.fprintf(file, "\n")


    for cell, i in mesh.cells {
    	fmt.fprintf(file, "%e, %e", cell.position.x, cell.position.y)
    	for vf in vfs {
    		v := cfd.vector_field_at(vf, i)
    		fmt.fprintf(file, ",%e,%e", v.x, v.y)
    	}
    	for sf in sfs {
    		fmt.fprintf(file, ",e", sf.data[i])
    	}
        fmt.fprintf(file, "\n")
    }
    return output_path, true
}