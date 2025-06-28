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
import "core:log"
import "core:os"

write_csv :: proc(mesh: cfd.Mesh, vfs: []Out_Vector_Field, sfs: []Out_Scalar_Field, path: string) -> (ok: bool) {
	file, err := os.open(path, os.O_CREATE | os.O_TRUNC | os.O_RDWR, 0o777)
	if err != nil {
		log.errorf("Could not create output CSV file %s", path)
		return false
	}
	defer os.close(file)


	fmt.fprintf(file, "x,y")


	for vf in vfs {
		fmt.fprintf(file, ",%s.x", vf.name)
		fmt.fprintf(file, ",%s.y", vf.name)
	}
	for sf in sfs {
		fmt.fprintf(file, ",%", sf.name)
	}
	fmt.fprintf(file, "\n")


	for cell, i in mesh.cells {
		fmt.fprintf(file, "%e, %e", cell.position.x, cell.position.y)
		for vf in vfs {
			fmt.fprintf(file, ",%e,%e", vf.data.x[i], vf.data.y[i])
		}
		for sf in sfs {
			fmt.fprintf(file, ",e", sf.data[i])
		}
		fmt.fprintf(file, "\n")
	}
	return true
}
