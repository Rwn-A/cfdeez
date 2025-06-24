package solvers

import "../../cfd"

import "core:log"

ITERS :: 10
MOMENTUM_RELAX :: 0.7
PRESSURE_RELAX :: 0.2

simple :: proc(
	mesh: cfd.Mesh,
	u: ^cfd.Vector_Field,
	p: ^cfd.Field,
	mass_flux: ^cfd.Face_Data, // pointer pass is not strictly nessecary, just to indicate it's modified.
	systems: [2]^cfd.Linear_System_Builder,
	viscosity, density: f64
) -> cfd.Linsolve_Error {
	solver_mem := cfd.solver_memory_create()
    defer cfd.solver_memory_destroy(&solver_mem)

    Ap := cfd.field_pool_cell_data(&cfd.field_pool)
	defer cfd.field_pool_return(&cfd.field_pool, Ap)

	u_old := cfd.field_pool_vector_cell_data(&cfd.field_pool)
	defer cfd.field_pool_return(&cfd.field_pool, ..u_old[:])

	H := cfd.field_pool_vector_cell_data(&cfd.field_pool)
	defer cfd.field_pool_return(&cfd.field_pool, ..H[:])

	p_old := cfd.field_pool_cell_data(&cfd.field_pool)
	defer cfd.field_pool_return(&cfd.field_pool, p_old)

	for iter in 0..<ITERS {
		cfd.linear_system_builder_reset(systems[0])
		cfd.linear_system_builder_reset(systems[1])
		for _, cell_id in mesh.cells {
            cfd.linear_system_builder_start_cell(systems[0], cell_id)
            cfd.linear_system_builder_start_cell(systems[1], cell_id)
            defer cfd.linear_system_builder_end_cell(systems[0])
            defer cfd.linear_system_builder_end_cell(systems[1])
            cfd.add_advection(systems[0], mesh, u.components.x, cfd.field_flux(mesh, u), cell_id)
            cfd.add_advection(systems[1], mesh, u.components.y, cfd.field_flux(mesh, u), cell_id)
            systems[0].mode = .Sub; systems[1].mode = .Sub
            cfd.add_laplacian(systems[0], mesh, &u.components.x, cell_id, viscosity)
            cfd.add_laplacian(systems[1], mesh, &u.components.y, cell_id, viscosity)
        }
		systems[0].mode = .Sub; systems[1].mode = .Sub
		grad := cfd.field_gradient(mesh, p)
		cfd.linear_system_builder_add_source(systems[0], mesh, cast([]f64)grad.x)
		cfd.linear_system_builder_add_source(systems[1], mesh, cast([]f64)grad.y)

		ls_x := cfd.linear_system_builder_assemble(systems[0])
	    ls_y := cfd.linear_system_builder_assemble(systems[1])

		cfd.linear_system_get_diagonal(ls_x, mesh, Ap)

		copy(u_old.x, u.components.x.data)
		copy(u_old.y, u.components.y.data)

		cfd.linear_system_solve_pgs(ls_x, &u.components.x, &solver_mem) or_return
	    cfd.linear_system_solve_pgs(ls_y, &u.components.y, &solver_mem) or_return

	    for _, cell_id in mesh.cells {
	    	xdat := &u.components.x.data[cell_id]
	    	ydat := &u.components.y.data[cell_id]
	    	xdat^ = MOMENTUM_RELAX * xdat^ + (1.0 - MOMENTUM_RELAX) * u_old.x[cell_id]
			ydat^ = MOMENTUM_RELAX * ydat^ + (1.0 - MOMENTUM_RELAX) * u_old.y[cell_id]
	    }

	    for _, cell_id in mesh.cells {
	        H.x[cell_id] = ls_x.b[cell_id]
	        H.y[cell_id] = ls_y.b[cell_id]
	        
	        for i in ls_x.M.row_indices[cell_id]..<ls_x.M.row_indices[cell_id + 1] {
	            if ls_x.M.columns[i] != cell_id {
	               H.x[cell_id] -= ls_x.M.values[i] * u.components.x.data[ls_x.M.columns[i]]
	            }
	        }
	        
	        for i in ls_y.M.row_indices[cell_id]..<ls_y.M.row_indices[cell_id + 1] {
	            if ls_y.M.columns[i] != cell_id {
	               	H.y[cell_id] -= ls_y.M.values[i] * u.components.y.data[ls_y.M.columns[i]]
	            }
	        }

	        H.x[cell_id] /= Ap[cell_id]  
	        H.y[cell_id] /= Ap[cell_id] 
	    }

	    // TODO: these BC's aren't completely accurate at outflow's because of non-zero pressure gradient
	    temp_H := cfd.Vector_Field{
	            components = {
	                {data = H.x, boundary_conditions = u.components.x.boundary_conditions},
	                {data = H.y, boundary_conditions = u.components.y.boundary_conditions},
	        },
	    }

	    cfd.linear_system_builder_reset(systems[0])
	    cfd.linear_system_builder_reset(systems[1])

	    copy(p_old, p.data)
	    for &p in p.data do p = 0

	    for _, cell_id in mesh.cells {
	        cfd.linear_system_builder_start_cell(systems[0], cell_id)
	        defer cfd.linear_system_builder_end_cell(systems[0])
	        cfd.add_laplacian(systems[0], mesh, p, cell_id, 1 / Ap[cell_id])
	    }
	    cfd.linear_system_builder_add_source(systems[0], mesh, cfd.field_divergence(mesh, &temp_H))

	    ls_p := cfd.linear_system_builder_assemble(systems[0])
	    
	    cfd.linear_system_solve_pcg(ls_p, p, &solver_mem) or_return
	    cfd.linear_system_builder_reset(systems[0])

	    grad = cfd.field_gradient(mesh, p)
	    for _, cell_id in mesh.cells {
	        d := density / Ap[cell_id]
	        u.components.x.data[cell_id] -= d * grad.x[cell_id]
	        u.components.y.data[cell_id] -= d * grad.y[cell_id]
	    }

	    avg_corr: f64
	    for dat in p.data {
	    	avg_corr += dat
	    }
	    avg_corr /= cast(f64)len(mesh.cells)

	    for _, cell_id in mesh.cells {
	        p.data[cell_id] = p_old[cell_id] + PRESSURE_RELAX * p.data[cell_id]
	    }

	    cfd.vfield_dirty(u)
	    copy(mass_flux^, cfd.field_flux(mesh, u))
	    for &flux in mass_flux do flux *= density

	   	if avg_corr < 1e-5 do break
	   	log.info(avg_corr)
	} 

		return .None
}
