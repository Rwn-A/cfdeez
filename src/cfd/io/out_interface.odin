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
import "core:path/filepath"

Out_Scalar_Field :: struct {
	data: []f64,
	name: string,
}

Out_Vector_Field :: struct {
	data: [cfd.VF_COMPONENTS][]f64,
	name: string,
}

Output_Fn :: #type proc(mesh: cfd.Mesh, vfs: []Out_Vector_Field, sfs: []Out_Scalar_Field, path: string) -> bool

output_path :: proc(dir: string, sim_name: string, ext: string, step: int) -> string {
	@(static) output_path_bufffer: [128]byte
	return fmt.bprintf(output_path_bufffer[:], "%s%c%s_%d.%s", dir, filepath.SEPARATOR, sim_name, step, ext)
}
