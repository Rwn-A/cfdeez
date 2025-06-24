package shared

import "core:os"
import "core:log"
import vmem"core:mem/virtual"
import sa"core:container/small_array"
import "core:path/filepath"
import "core:math"
import "core:strings"
import "core:fmt"
import "core:slice"

import "../fml"
import "../cfd"
import cfd_io"../cfd/io"

// TODO: validate config options make sense.
//TODO: validate boundaries

Physics_Option :: enum {Transport, IncFlow}
Output_Option :: enum {CSV, VTU}

// This is the head honcho of the simulation.
// It has all the configuration and memory needed to execute the simulation.
Case :: struct{
    arena: ^vmem.Arena,
    name: string,
    mesh: cfd.Mesh,
    fields: cfd.Primary_Fields,
    fluid: struct {density, viscosity: f64},
    physics: bit_set[Physics_Option],
    output: struct {
        directory: string,
        format: []cfd_io.Output_Fn,
        pvd: bool,
        frequency: int,
    },
    timestep: Maybe(f64),
    steps: int,
    passive_names: []string, //used for output files, index here matches index in fields.passives array.
    passive_diffusivity: []f64,


}

//fixed value or expr
Inflow_Profile_Schema :: Maybe(union{f64, fml.Expression})

//fixed value, file or expr
Initial_Condition_Schema :: Maybe(union{f64, string, fml.Expression})

Passive_Schema :: struct {
    name: string,
    diffusivity: f64,
    inflow_profile: Inflow_Profile_Schema,
    initial_condition: Initial_Condition_Schema,
}

Case_Schema :: struct {
    name: string,
    mesh: struct {path: string},
    fluid: struct {density, viscosity: f64},
    physics: []Physics_Option,
    output: struct {
        directory: string,
        formats: []Output_Option,
    },
    time: Maybe(struct {
        timestep: f64,
        steps: int,
        output_frequency: Maybe(int),
        enable_pvd: Maybe(bool),
    }),
    boundaries: struct {
        wall: Maybe([]string),
        inflow: Maybe([]string),
        outflow: Maybe([]string),
    },
    velocity: struct {
       inflow_profile: [2]Inflow_Profile_Schema,
       initial_conditions: [2]Initial_Condition_Schema,
    },
    passives: Maybe([]Passive_Schema),
}

load_case :: proc(path: string, arr: ^vmem.Arena) -> (c: Case, success: bool) {
    scratch := vmem.Arena{}
    if err := vmem.arena_init_growing(&scratch); err != nil {
        log.fatal("Unable to allocate")
        return {}, false
    }
    defer vmem.arena_destroy(&scratch)
    scratch_allocator := vmem.arena_allocator(&scratch)

    content, ok := os.read_entire_file_from_filename(path, scratch_allocator)
    if !ok {
        log.errorf("File %s does not exist or can't be read.", path)
        return {}, false
    }

    schema: Case_Schema
    if err := fml.unmarshal(content, path, &schema, scratch_allocator); err != nil {
        fml.log_error(err)
        return {}, false
    }

    os.set_current_directory(filepath.dir(path))

    c.arena = arr
    if err := vmem.arena_init_growing(c.arena); err != nil {
        log.fatal("Unable to allocate.")
        return {}, false
    }
    defer if success == false do vmem.arena_destroy(c.arena)

    context.allocator = vmem.arena_allocator(c.arena)

    c.name = strings.clone(schema.name)
    c.fluid = schema.fluid
    for opt in schema.physics do c.physics += {opt}
    for opt in schema.output.formats {
        builder := make([dynamic]cfd_io.Output_Fn)
        switch opt {
            case .CSV: append(&builder, cfd_io.to_csv)
            case .VTU:  append(&builder, cfd_io.to_vtu)
        }
        c.output.format = builder[:]
    }
    c.output.directory = strings.clone(schema.output.directory)
    if time, has_time := schema.time.?; has_time {
        c.timestep = time.timestep
        c.steps = time.steps
        c.output.frequency = time.output_frequency.? or_else 1
        c.output.pvd = time.enable_pvd.? or_else true
    }else{
        c.steps = 1
        c.output.frequency = 1
    }

    if mesh_fd, err := os.open(schema.mesh.path); err == nil {
        defer os.close(mesh_fd)
        c.mesh = cfd_io.from_gmsh(os.stream_from_handle(mesh_fd), c.arena) or_return
    }else{
        log.errorf("Failed to open mesh file %s, does it exist?", schema.mesh.path)
        return {}, false
    }

    cfd.field_pool_init(&cfd.field_pool, c.mesh)

    empty_field :: proc() -> cfd.Field {
        return {data = cfd.field_pool_cell_data(&cfd.field_pool), boundary_conditions = make(cfd.Boundaries)}
    }

    c.fields.u.components = {
        empty_field(),
        empty_field(),
    }
    c.fields.p = empty_field()

    load_ic(c.mesh, c.fields.u.components.x, schema.velocity.initial_conditions.x)
    load_ic(c.mesh, c.fields.u.components.y,  schema.velocity.initial_conditions.y)

    field_builder := make([dynamic]cfd.Field)
    name_builder := make([dynamic]string)
    diffusivity_builder := make([dynamic]f64)
    if passives, has := schema.passives.?; has {
        for passive in passives {
            field := empty_field()
            append(&name_builder, strings.clone(passive.name))
            append(&diffusivity_builder, passive.diffusivity)
            load_ic(c.mesh, field, passive.initial_condition) or_return
            append(&field_builder, field)
        }
        c.fields.passives = field_builder[:]
        c.passive_names = name_builder[:]
        c.passive_diffusivity = diffusivity_builder[:]
    }


    if wall_bnds, has := schema.boundaries.wall.?; has {
        for bnd_name in wall_bnds {
            cfd.field_set_bnd_from_name(&c.fields.u.components.x, c.mesh, bnd_name, {.Dirichlet, 0})
            cfd.field_set_bnd_from_name(&c.fields.u.components.y, c.mesh, bnd_name, {.Dirichlet, 0})
            cfd.field_set_bnd_from_name(&c.fields.p, c.mesh, bnd_name, {.Neumann, 0})
            for &passive in c.fields.passives {
                cfd.field_set_bnd_from_name(&passive, c.mesh, bnd_name, {.Dirichlet, 0})
            }
        }
    }

    if outflow_bnd, has := schema.boundaries.outflow.?; has {
        for bnd_name in outflow_bnd {
            cfd.field_set_bnd_from_name(&c.fields.u.components.x, c.mesh, bnd_name, {.Neumann, 0})
            cfd.field_set_bnd_from_name(&c.fields.u.components.y, c.mesh, bnd_name, {.Neumann, 0})
            cfd.field_set_bnd_from_name(&c.fields.p, c.mesh, bnd_name, {.Dirichlet, 0})
            for &passive in c.fields.passives {
                cfd.field_set_bnd_from_name(&passive, c.mesh, bnd_name, {.Neumann, 0})
            }
        }
    }

    if  inflow_bnd, has := schema.boundaries.inflow.?; has {
        for bnd_name in inflow_bnd {
            //split up the inflow boundary into pieces so we can assign the profile to it.
            //this is scuffed
            //this also invalidates the mesh.boundary_names, but maybe thats what we want so it still matches the input mesh.
            first_new_bnd_id: int
            original_inflow_id: int
            for k, v in c.mesh.boundary_names {
                if k == bnd_name do original_inflow_id = v
                if v > first_new_bnd_id do first_new_bnd_id = v
            }
            first_new_bnd_id += 1

            for &face, i in c.mesh.faces {
                if face.boundary_tag != original_inflow_id do continue
                new_id := first_new_bnd_id + i
                face.boundary_tag = first_new_bnd_id + i
                cfd.field_set_bnd(&c.fields.u.components.x, new_id, {.Dirichlet, inflow_value(schema.velocity.inflow_profile.x, face) or_return})
                cfd.field_set_bnd(&c.fields.u.components.y, new_id, {.Dirichlet, inflow_value(schema.velocity.inflow_profile.y, face) or_return})
                cfd.field_set_bnd(&c.fields.p, new_id, {.Neumann, 0})
                for &passive, i in c.fields.passives {
                    value := inflow_value(schema.passives.?[i].inflow_profile, face) or_return
                    cfd.field_set_bnd(&passive, new_id, {.Dirichlet, value})
                }
            }
        }
    }

    // Validation

    if c.output.pvd {
        if _, has_timestep := c.timestep.?; !has_timestep {
            log.error("PVD files can only be generated for transient simulations.")
            return {}, false
        }
        if !slice.contains(schema.output.formats, Output_Option.VTU) {
            log.error("PVD Files can only be generated from VTU files, include VTU in your output formats.")
            return {}, false
        }
    }


    return c, true

    load_ic :: proc(mesh: cfd.Mesh, field: cfd.Field, ic_schema: Initial_Condition_Schema) -> bool {
        if _, has := ic_schema.?; !has do return true
        switch value in ic_schema.? {
            case string: unimplemented("Loading initial conditions from files.")
            case f64:
                for &d in field.data do d = value
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
            case f64: return value, true
            case fml.Expression:
                host_vars := make(map[string]f64, context.temp_allocator)
                defer free_all(context.temp_allocator)
                host_vars["x"] = face.position.x
                host_vars["y"] = face.position.y
                return fml.evaluate(value, host_vars, expr_host_fn)
            case: unreachable()
        }

    }
}

expr_host_fn :: proc(name: string, stack: ^fml.Eval_Stack) -> bool {
    switch name {
        case "cos": sa.push_back(stack, math.cos(sa.pop_back(stack)))
        case "sin": sa.push_back(stack, math.sin(sa.pop_back(stack)))
        case "tan": sa.push_back(stack, math.tan(sa.pop_back(stack)))
        case "exp": sa.push_back(stack, math.exp(sa.pop_back(stack)))
        case "pow": sa.push_back(stack, math.pow(sa.pop_back(stack), sa.pop_back(stack)))
        case "ln":  sa.push_back(stack, math.ln(sa.pop_back(stack)))
        case "log": sa.push_back(stack, math.log(sa.pop_back(stack), sa.pop_back(stack)))
        case:
            log.errorf("function %s is unknown.", name); return false
    }
    return true
}
