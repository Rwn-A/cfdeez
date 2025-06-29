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
package app


import sa "core:container/small_array"
import "core:fmt"
import "core:log"
import "core:math"
import vmem "core:mem/virtual"
import "core:os"
import "core:path/filepath"
import "core:slice"
import "core:strings"

import "../cfd"
import cfd_io "../cfd/io"
import "../fml"

Case :: struct {
	arena:                 vmem.Arena,
	fp:                    cfd.Field_Data_Pool,
	name:                  string,
	mesh:                  cfd.Mesh,
	fields:                Primary_Fields,
	fluid:                 struct {
		density, viscosity: f64,
	},
	physics:               bit_set[Physics_Option],
	output:                struct {
		directory: string,
		format:   bit_set[Output_Option],
		frequency: int,
	},
	timestep:              Maybe(f64),
	steps:                 int,
	passive_names:         []string, //used for output files, index here matches index in fields.passives array.
	passive_diffusivities: []f64, //index here matches index in fields.passives array.
}

Primary_Fields :: struct {
	u:        cfd.Vector_Field,
	p:        cfd.Scalar_Field,
	passives: []cfd.Scalar_Field,
}

Physics_Option :: enum {
	Transport,
	IncFlow,
}

Output_Option :: enum {
	CSV,
	VTK,
}

//fixed value or expr
Inflow_Profile_Schema :: Maybe(union {
		f64,
		fml.Expression,
	})

//fixed value, file or expr
Initial_Condition_Schema :: Maybe(union {
		f64,
		string,
		fml.Expression,
	})



Passive_Schema :: struct {
	name:              string,
	diffusivity:       f64,
	inflow_profile:    Inflow_Profile_Schema,
	initial_condition: Initial_Condition_Schema,
}

Case_Schema :: struct {
	name:       string,
	mesh:       struct {
		path: string,
	},
	fluid:      struct {
		density, viscosity: f64,
	},
	physics:    []Physics_Option,
	output:     struct {
		directory: string,
		formats:   []Output_Option,
	},
	time:       Maybe(struct {
			timestep:         f64,
			steps:            int,
			output_frequency: Maybe(int),
		}),
	boundaries: struct {
		wall:    Maybe([]string),
		inflow:  Maybe([]string),
		outflow: Maybe([]string),
	},
	velocity:   struct {
		inflow_profile:     [2]Inflow_Profile_Schema,
		initial_conditions: [2]Initial_Condition_Schema,
	},
	passives:   Maybe([]Passive_Schema),
}

load_case :: proc(c: ^Case, path: string) -> (success: bool) {
	scratch := vmem.Arena{}
	if err := vmem.arena_init_growing(&scratch); err != nil {
		log.fatal("Unable to allocate")
		return false
	}
	defer vmem.arena_destroy(&scratch)
	scratch_allocator := vmem.arena_allocator(&scratch)

	content, ok := os.read_entire_file_from_filename(path, scratch_allocator)
	if !ok {
		log.errorf("File %s does not exist or can't be read.", path)
		return false
	}

	schema: Case_Schema
	if err := fml.unmarshal(content, path, &schema, scratch_allocator); err != nil {
		fml.log_error(err)
		return false
	}

	config_directory := filepath.dir(path, scratch_allocator)
	os.set_current_directory(config_directory)

	if err := vmem.arena_init_growing(&c.arena); err != nil {
		log.fatal("Unable to allocate.")
		return false
	}
	defer if success == false do vmem.arena_destroy(&c.arena)

	context.allocator = vmem.arena_allocator(&c.arena)

	c.name = strings.clone(schema.name)

	c.fluid = schema.fluid

	for opt in schema.physics {c.physics += {opt}}
	for opt in schema.output.formats {c.output.format += {opt}}

	c.output.directory = strings.clone(schema.output.directory)

	if time, has_time := schema.time.?; has_time {
		c.timestep = time.timestep
		c.steps = time.steps
		c.output.frequency = time.output_frequency.? or_else 1
	} else {
		c.steps = 1
		c.output.frequency = 1
	}

	bnd_map: map[string]cfd.BoundaryId
	if mesh_fd, err := os.open(schema.mesh.path); err == nil {
		defer os.close(mesh_fd)
		c.mesh, bnd_map = cfd_io.from_gmsh(os.stream_from_handle(mesh_fd), &c.arena) or_return
	} else {
		log.errorf("Failed to open mesh file %s, does it exist?", schema.mesh.path)
		return false
	}

	cfd.field_pool_init(&c.fp, c.mesh)

	c.fields.u = cfd.vector_field_new(&c.fp)
	c.fields.p = cfd.scalar_field_new(&c.fp)
	cfd.field_start_mutation(&c.fields.u)
	cfd.field_start_mutation(&c.fields.p)
	defer cfd.field_end_mutation(&c.fields.u)
	defer cfd.field_end_mutation(&c.fields.p)

	load_ic(c.mesh, c.fields.u.components.x, schema.velocity.initial_conditions.x)
	load_ic(c.mesh, c.fields.u.components.y, schema.velocity.initial_conditions.y)

	//load passives
	field_builder := make([dynamic]cfd.Scalar_Field)
	name_builder := make([dynamic]string)
	diffusivity_builder := make([dynamic]f64)
	if passives, has := schema.passives.?; has {
		for passive in passives {
			field := cfd.scalar_field_new(&c.fp)
			cfd.field_start_mutation(&field)
			append(&name_builder, strings.clone(passive.name))
			append(&diffusivity_builder, passive.diffusivity)
			load_ic(c.mesh, field, passive.initial_condition) or_return
			cfd.field_end_mutation(&field)
			append(&field_builder, field)
		}
		c.fields.passives = field_builder[:]
		c.passive_names = name_builder[:]
		c.passive_diffusivities = diffusivity_builder[:]
	}

	//boundary conditions
	if wall_bnds, has := schema.boundaries.wall.?; has {
		for bnd_name in wall_bnds {
			c.fields.u.components.x.bnds[bnd_map[bnd_name]] = {.Dirichlet, 0}
			c.fields.u.components.y.bnds[bnd_map[bnd_name]] = {.Dirichlet, 0}
			c.fields.p.bnds[bnd_map[bnd_name]] = {.Neumann, 0}
			for &passive in c.fields.passives {
				passive.bnds[bnd_map[bnd_name]] = {.Dirichlet, 0}
			}
		}
	}


	if outflow_bnds, has := schema.boundaries.outflow.?; has {
		for bnd_name in outflow_bnds {
			c.fields.u.components.x.bnds[bnd_map[bnd_name]] = {.Neumann, 0}
			c.fields.u.components.y.bnds[bnd_map[bnd_name]] = {.Neumann, 0}
			c.fields.p.bnds[bnd_map[bnd_name]] = {.Dirichlet, 0}
			for &passive in c.fields.passives {
				passive.bnds[bnd_map[bnd_name]] = {.Neumann, 0}
			}
		}
	}

	if inflow_bnd, has := schema.boundaries.inflow.?; has {
		for bnd_name in inflow_bnd {
			//split up the inflow boundary into pieces so we can assign the profile to it.
			//this is scuffed
			first_new_bnd_id: cfd.BoundaryId
			original_inflow_id: cfd.BoundaryId
			for k, v in bnd_map {
				if k == bnd_name do original_inflow_id = v
				if v > first_new_bnd_id do first_new_bnd_id = v
			}
			first_new_bnd_id += 1

			for &face, i in c.mesh.faces {
				if id, has := face.secondary.(cfd.BoundaryId); has {
					if id != original_inflow_id {continue}
				}else{
					continue
				}
				new_id := first_new_bnd_id + cfd.BoundaryId(i)
				face.secondary = new_id
				c.fields.u.components.x.bnds[new_id] = {
					.Dirichlet,
					inflow_value(schema.velocity.inflow_profile.x, face) or_return,
				}
				c.fields.u.components.y.bnds[new_id] = {
					.Dirichlet,
					inflow_value(schema.velocity.inflow_profile.y, face) or_return,
				}
				c.fields.p.bnds[new_id] = {.Neumann, 0}

				for &passive, i in c.fields.passives {
					passive.bnds[new_id] = {.Dirichlet, inflow_value(schema.passives.?[i].inflow_profile, face) or_return}
				}
			}
		}
	}

	load_ic :: proc(mesh: cfd.Mesh, field: cfd.Scalar_Field, ic_schema: Initial_Condition_Schema) -> bool {
		if _, has := ic_schema.?; !has do return true
		switch value in ic_schema.? {
		case string:
			unimplemented("Loading initial conditions from files.")
		case f64:
			for &d in field.data {d = value}
		case fml.Expression:
			host_vars := make(map[string]f64, context.temp_allocator)
			defer free_all(context.temp_allocator)
			for cell, idx in mesh.cells {
				host_vars["x"] = cell.position.x
				host_vars["y"] = cell.position.y
				field.data[idx] = fml.evaluate(value, host_vars, expr_host_fn) or_return
			}
		}
		return true
	}

	inflow_value :: proc(schema: Inflow_Profile_Schema, face: cfd.Face) -> (f: f64, ok: bool) {
		if _, has := schema.?; !has {
			log.error("Inflow boundary specified in boundaries table, but inflow profile was not supplied for every field.")
			return 0, false
		}

		switch value in schema.? {
		case f64:
			return value, true
		case fml.Expression:
			host_vars := make(map[string]f64, context.temp_allocator)
			defer free_all(context.temp_allocator)
			host_vars["x"] = face.position.x
			host_vars["y"] = face.position.y
			return fml.evaluate(value, host_vars, expr_host_fn)
		case:
			unreachable()
		}

	}

	return true
}


expr_host_fn :: proc(name: string, stack: ^fml.Eval_Stack) -> bool {
	switch name {
	case "cos":
		sa.push_back(stack, math.cos(sa.pop_back(stack)))
	case "sin":
		sa.push_back(stack, math.sin(sa.pop_back(stack)))
	case "tan":
		sa.push_back(stack, math.tan(sa.pop_back(stack)))
	case "exp":
		sa.push_back(stack, math.exp(sa.pop_back(stack)))
	case "pow":
		sa.push_back(stack, math.pow(sa.pop_back(stack), sa.pop_back(stack)))
	case "ln":
		sa.push_back(stack, math.ln(sa.pop_back(stack)))
	case "log":
		sa.push_back(stack, math.log(sa.pop_back(stack), sa.pop_back(stack)))
	case:
		log.errorf("function %s is unknown.", name);return false
	}
	return true
}
