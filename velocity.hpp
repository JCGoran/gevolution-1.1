#include "psi.h"
#include "geometry.h"
#include "refine.h"
#include "grid.h"

psi_int psi_aabb_periodic(psi_rvec* pos, psi_rvec* rbox, psi_rvec* window, psi_rvec* box);

template<typename part, typename part_info, typename part_dataType>
void test_id(Particles<part, part_info, part_dataType> * pcls, long *ID, double boxsize, Field<Real> *density, Field<Real> *velocity)
{
    Site xPart(pcls->lattice());
    typename std::list<part>::iterator it;

    psi_grid grid;
    memset(&grid, 0, sizeof(grid));
    grid.type = 0;
    grid.dim = 3;
    grid.fields[0] = density->data();
    grid.fields[1] = velocity->data();
    grid.n.i = pcls->lattice().localSize(0);
    grid.n.j = pcls->lattice().localSize(1);
    grid.n.k = pcls->lattice().localSize(2);
    grid.window[0].x = 0.;
    grid.window[0].y = 1.0*coordSkip(0)/size(1);
    grid.window[0].z = 1.0*coordSkip(1)/size(2);
    grid.window[1].x = 1.;
    grid.window[1].y = 1.0*(coordSkip(0) + sizeLocal(1))/size(1);
    grid.window[1].z = 1.0*(coordSkip(1) + sizeLocal(2))/size(2);
    grid.d.x = (grid.window[1].x - grid.window[0].x)/grid.n.i;
    grid.d.y = (grid.window[1].y - grid.window[0].y)/grid.n.j;
    grid.d.z = (grid.window[1].z - grid.window[0].z)/grid.n.k;

	psi_int e, nghosts, tind, internal_rtree, e1, t1;

	// a local copy of the position in case it is modified due to periodicity
	psi_int vpere = 8;
	psi_rvec tpos[8], tvel[8], tpos1[8], tvel1[8];
	psi_real tmass, tmass1;
	psi_rvec trbox[2];
	psi_rvec *pos0, *vel0, *pos1, *vel1;
	psi_rvec rbox0[2], rbox1[2], mbox[2];
	psi_real mass0, mass1;

    const int meshdim = 3;

	// ghost tetrahedra for handling periodicity
	psi_rvec gpos[32];
	psi_rvec grbox[8];

	// constant-size buffer for max refinement level
	psi_rtree_query qry;
	psi_tet_buffer tetbuf, tetbuf1;
	psi_tet_buffer_init(&tetbuf, 0, 0);

    psi_rvec box[2];
    box[0].x = 0.;
    box[0].y = 0.;
    box[0].z = 0.;
    box[1].x = 1.;
    box[1].y = 1.;
    box[1].z = 1.;



    psi_rvec pos[8], vel[8];

    for (xPart.first(); xPart.test(); xPart.next()){
        for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it){
            for (psi_int i = 0; i < 8; ++i){
                // find the right 8 particles
                if (it->ID == ID[i]){
                    pos[i].x = it->pos[0];
                    pos[i].y = it->pos[1];
                    pos[i].z = it->pos[2];
                    vel[i].x = it->vel[0];
                    vel[i].y = it->vel[1];
                    vel[i].z = it->vel[2];
                }
            }
        }
    }
    for (xPart.first(); xPart.test(); xPart.next()){
        for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it){

            // a local copy of the element, in case it's modified
            // make it periodic, get its bounding box, and check it against the grid
            tmass = 1;//mesh->mass[e]; 
            if(!psi_aabb_periodic(pos, trbox, grid.window, box)){
                continue;
            }

            // refine elements into the tet buffer 
            psi_tet_buffer_refine(&tetbuf, pos, vel, tmass, 8);

            // loop over each tet in the buffer
            // copy the position to the ghost array and compute its aabb
            // make ghosts and sample tets to the grid
            for(psi_int t = 0; t < tetbuf.num; ++t){
                memcpy(gpos, &tetbuf.pos[(meshdim+1)*t], (meshdim+1)*sizeof(psi_rvec));
                psi_aabb(gpos, meshdim+1,grbox);
                //psi_make_ghosts(gpos, grbox, &nghosts, (meshdim+1), grid->window, mesh);
                nghosts = 1;
                for(psi_int g = 0; g < nghosts; ++g){
                    psi_voxelize_tet(&gpos[(meshdim+1)*g], &tetbuf.vel[(meshdim+1)*t], tetbuf.mass[t], &grbox[2*g], &grid);
                }
            }
        }
    }
}


psi_int psi_aabb_periodic(psi_rvec* pos, psi_rvec* rbox, psi_rvec* window, psi_rvec* box) {

	// pos and rbox can both be modified!
	psi_int i, v;
	psi_real span;

	psi_aabb(pos, 8, rbox);

	// if the grid is periodic, move the vertices as needed
    for(i = 0; i < 3; ++i) {
        span = box[1].xyz[i] - box[0].xyz[i];
        if(rbox[1].xyz[i] - rbox[0].xyz[i] > 0.5*span) {
            rbox[0].xyz[i] = 1.0e30;
            rbox[1].xyz[i] = -1.0e30;
            for(v = 0; v < 8; ++v) {
                // wrap elements that cross the right-hand boundary while recomputing the box
                if(pos[v].xyz[i] > box[0].xyz[i] + 0.5*span) pos[v].xyz[i] -= span;
                if(pos[v].xyz[i] < rbox[0].xyz[i]) rbox[0].xyz[i] = pos[v].xyz[i];
                if(pos[v].xyz[i] > rbox[1].xyz[i]) rbox[1].xyz[i] = pos[v].xyz[i];
            }
        }
        if(rbox[0].xyz[i] > window[1].xyz[i]) return 0; 
        if(rbox[1].xyz[i] < window[0].xyz[i]
            && rbox[0].xyz[i]+span > window[1].xyz[i]) return 0; 
    }
	return 1;
}

/***
void psi_voxels(psi_grid* grid, Particles<part, part_info, part_dataTyle> *pcls)
{
	//setbuf(stdout, NULL);
	psi_int e, g, t, v, nghosts, tind, internal_rtree, e1, t1;

	// a local copy of the position in case it is modified due to periodicity
	psi_int vpere = mesh->elemtype;
	psi_rvec tpos[vpere], tvel[vpere], tpos1[vpere], tvel1[vpere];
	psi_real tmass, tmass1;
	psi_rvec trbox[2];
	psi_rvec *pos0, *vel0, *pos1, *vel1;
	psi_rvec rbox0[2], rbox1[2], mbox[2];
	psi_real mass0, mass1;

	// ghost tetrahedra for handling periodicity
	psi_rvec gpos[(1<<mesh->dim)*(mesh->dim+1)];
	psi_rvec grbox[(1<<mesh->dim)*2];

	// constant-size buffer for max refinement level
	psi_rtree_query qry;
	psi_tet_buffer tetbuf, tetbuf1;
	psi_tet_buffer_init(&tetbuf, reftol, max_ref_lvl);

	// check to see if we must build an rtree for this function
	// if so, create one from the mesh


	// loop over all elements
    Site xPart(pcls->lattice());
    typename std::list<part>::iterator it;
    for (xPart.first(); xPart.test(); xPart.next()){
        for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it){
            //printf("part ID: %ld\n", it->ID);
        }
    }

		// a local copy of the element, in case it's modified
		// make it periodic, get its bounding box, and check it against the grid
		tmass = mesh->mass[e]; 
		for(v = 0; v < vpere; ++v) {
			tind = mesh->connectivity[e*vpere+v];
			tpos[v] = mesh->pos[tind];
			tvel[v] = mesh->vel[tind];
		}
		if(!psi_aabb_periodic(tpos, trbox, grid->window, mesh)) 
			continue;

		// refine elements into the tet buffer 
		psi_tet_buffer_refine(&tetbuf, tpos, tvel, tmass, mesh->elemtype);

		// linear in the density
		if(mode == PSI_MODE_DENSITY) {

			// loop over each tet in the buffer
			// copy the position to the ghost array and compute its aabb
			// make ghosts and sample tets to the grid
			for(t = 0; t < tetbuf.num; ++t) {
				memcpy(gpos, &tetbuf.pos[(mesh->dim+1)*t], (mesh->dim+1)*sizeof(psi_rvec));
				psi_aabb(gpos, mesh->dim+1,grbox);
				psi_make_ghosts(gpos, grbox, &nghosts, (mesh->dim+1), grid->window, mesh);
				for(g = 0; g < nghosts; ++g) 
					psi_voxelize_tet(&gpos[(mesh->dim+1)*g], &tetbuf.vel[(mesh->dim+1)*t], tetbuf.mass[t], &grbox[2*g], grid);
			}
		}

		// quadratic in the density
		else if(mode == PSI_MODE_ANNIHILATION) {

			// TODO: make and save the clip planes, density, rboxes for tetbuf 1 here

			// query the rtree for overlapping elements,
			// refine them, and intersect all pairs of resulting tets
			// TODO: ghosts and periodicity
			psi_rtree_query_init(&qry, rtree, trbox);
			while(psi_rtree_query_next(&qry, &e1)) {

				// get the overlapping element and refine it
				tmass1 = mesh->mass[e1]; 
				for(v = 0; v < vpere; ++v) {
					tind = mesh->connectivity[e1*vpere+v];
					tpos1[v] = mesh->pos[tind];
					tvel1[v] = mesh->vel[tind];
				}
				psi_tet_buffer_refine(&tetbuf1, tpos1, tvel1, tmass1, mesh->elemtype);

				// loop over each pair of tets between the refine buffers 
				for(t = 0; t < tetbuf.num; ++t) {

					// get the bounding box of the first tet
					// to use for the rtree query to find overlapping tets
					pos0 = &tetbuf.pos[(PSI_NDIM+1)*t];
					vel0 = &tetbuf.vel[(PSI_NDIM+1)*t];
					mass0 = tetbuf.mass[t];
					psi_aabb(pos0, PSI_NDIM+1, rbox0); 

					// loop over the second buffer
					for(t1 = 0; t1 < tetbuf1.num; ++t1) {
						pos1 = &tetbuf1.pos[(PSI_NDIM+1)*t1];
						vel1 = &tetbuf1.vel[(PSI_NDIM+1)*t1];
						mass1 = tetbuf1.mass[t1];
						psi_aabb(pos1, PSI_NDIM+1, rbox1); 

						// voxelize the annihilation rate between pairs of tets
						// only of mbox (the mutual intersection of bounding boxes) exists
						if(psi_aabb_ixn(rbox0, rbox1, mbox))
							psi_voxelize_annihilation(pos0, vel0, mass0, pos1, vel1, mass1, mbox, grid);
					}
				}
			}
		}

#ifdef PYMODULE
		if(PyErr_CheckSignals() < 0)
			return;
#endif
	}
	psi_printf("\n");

	// free temporary buffers
	psi_tet_buffer_destroy(&tetbuf);
	if(mode == PSI_MODE_ANNIHILATION)
		psi_tet_buffer_destroy(&tetbuf1);
	if(internal_rtree) 
		psi_rtree_destroy(rtree);

}

template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::coutPart(long ID)
{
    Site x(lat_part_);
    typename std::list<part>::iterator it;

    for(x.first();x.test();x.next())
    {
        if((field_part_)(x).size!=0)
        {
            for(it=(field_part_)(x).parts.begin(); it != (field_part_)(x).parts.end();++it)
            {
                if((*it).ID==ID)cout<< "Parallel ranks: ("<<parallel.grid_rank()[1]<<","<<parallel.grid_rank()[0]<<") ; "<< part_global_info_.type_name<<": "<<*it<< " , lattice position:"<<x <<endl;
            }
        }

    }

}


void psi_make_ghosts(psi_rvec* elems, psi_rvec* rboxes, psi_int* num, psi_int stride, psi_rvec* window, psi_mesh* mesh) {

	psi_int i, p, v, pflags;
	psi_real pshift;
	psi_rvec span;

	// if periodic, make ghosts on the fly
	(*num) = 1;
	pflags = 0;
	if(mesh->periodic) {
		for(i = 0; i < mesh->dim; ++i) {
			span.xyz[i] = mesh->box[1].xyz[i]-mesh->box[0].xyz[i]; 
			if(rboxes[0].xyz[i] < mesh->box[0].xyz[i]) pflags |= (1 << i);
		}
		for(p = 1; p < (1<<mesh->dim); ++p) {
			if((pflags & p) == p) {
				for(i = 0; i < mesh->dim; ++i) {
					pshift = span.xyz[i]*(1 & (p >> i));
					rboxes[2*(*num)+0].xyz[i] = rboxes[0].xyz[i] + pshift; 
					rboxes[2*(*num)+1].xyz[i] = rboxes[1].xyz[i] + pshift;
					if(rboxes[2*(*num)+0].xyz[i] > window[1].xyz[i]) goto next_shift;
					if(rboxes[2*(*num)+1].xyz[i] < window[0].xyz[i]) goto next_shift;
					for(v = 0; v < stride; ++v)
						elems[stride*(*num)+v].xyz[i] = elems[v].xyz[i] + pshift; 
				}
				(*num)++;
			}
			next_shift: continue;
		}	
	}
}

void psi_voxelize_tet(psi_rvec* pos, psi_rvec* vel, psi_real mass, psi_rvec* rbox, psi_rvec *window, int *griddim) {

	psi_int i, j, ii, jj, gridind, polyorder;
	psi_real cm[PSI_NDIM], cv[PSI_NDIM], sxx[PSI_NDIM][PSI_NDIM], sxv[PSI_NDIM][PSI_NDIM], svv[PSI_NDIM][PSI_NDIM];
	psi_real cwght, dwght;
	psi_real DX[PSI_NDIM][PSI_NDIM], DV[PSI_NDIM][PSI_NDIM];
	psi_real moments[10];
	psi_poly curpoly;

#if PSI_NDIM == 3
	const static psi_int qinds[3][3] = {{4,5,6},{5,7,8},{6,8,9}};
#elif PSI_NDIM == 2
	const static psi_int qinds[2][2] = {{3,4},{4,5}};
#endif

	// deformation matrices for the input tet's mass coordinates
	for(i = 0; i < PSI_NDIM; ++i)
	for(j = 0; j < PSI_NDIM; ++j) {
		DX[i][j] = pos[j+1].xyz[i]-pos[0].xyz[i];
		DV[i][j] = vel[j+1].xyz[i]-vel[0].xyz[i];
	}

	// initialize the tet as an edge-vertex graph 
	// and clamp it to the grid
	psi_init_tet(&curpoly, pos);

	// determine the correct poly order
	polyorder = 0;
	if(grid->fields[PSI_GRID_X] || grid->fields[PSI_GRID_V]) {
		polyorder = 1;	
	}
	if(grid->fields[PSI_GRID_XX] || grid->fields[PSI_GRID_XV] || grid->fields[PSI_GRID_VV]) {
		polyorder = 2;	
	}

	// initialize the voxelization iterator
	psi_voxels vox;
	psi_voxels_init(&vox, &curpoly, polyorder, rbox, grid);
	while(psi_voxels_next(&vox, moments, &gridind)) {
	
		// get moments for this fragment from the mass coordinates 
		if(grid->fields[PSI_GRID_X])
			for(i = 0; i < PSI_NDIM; ++i) {
				cm[i] = moments[0]*pos[0].xyz[i];
				for(j = 0; j < PSI_NDIM; ++j)
					cm[i] += moments[j+1]*DX[i][j];
				cm[i] /= moments[0];
			}
		if(grid->fields[PSI_GRID_V])
			for(i = 0; i < PSI_NDIM; ++i) {
				cv[i] = moments[0]*vel[0].xyz[i];
				for(j = 0; j < PSI_NDIM; ++j) 
					cv[i] += moments[j+1]*DV[i][j];
				cv[i] /= moments[0];
			}
		if(grid->fields[PSI_GRID_XX]) {
			for(i = 0; i < PSI_NDIM; ++i) 
			for(j = 0; j < PSI_NDIM; ++j)  {
				sxx[i][j] = moments[0]*(pos[0].xyz[i]-cm[i])*(pos[0].xyz[j]-cm[j]);
				for(ii = 0; ii < PSI_NDIM; ++ii) {
					sxx[i][j] += DX[i][ii]*moments[1+ii]*(pos[0].xyz[j]-cm[j]);
					sxx[i][j] += DX[j][ii]*moments[1+ii]*(pos[0].xyz[i]-cm[i]);
				}
				for(ii = 0; ii < PSI_NDIM; ++ii)
				for(jj = 0; jj < PSI_NDIM; ++jj)
					sxx[i][j] += DX[i][ii]*DX[j][jj]*moments[qinds[ii][jj]];
				sxx[i][j] /= moments[0];
			}
		}
		if(grid->fields[PSI_GRID_XV]) {
			for(i = 0; i < PSI_NDIM; ++i) 
			for(j = 0; j < PSI_NDIM; ++j)  {
				sxv[i][j] = moments[0]*(pos[0].xyz[i]-cm[i])*(vel[0].xyz[j]-cv[j]);
				for(ii = 0; ii < PSI_NDIM; ++ii) {
					sxv[i][j] += DX[i][ii]*moments[1+ii]*(pos[0].xyz[j]-cm[j]);
					sxv[i][j] += DV[j][ii]*moments[1+ii]*(vel[0].xyz[i]-cv[i]);
				}
				for(ii = 0; ii < PSI_NDIM; ++ii)
				for(jj = 0; jj < PSI_NDIM; ++jj)
					sxv[i][j] += DX[i][ii]*DV[j][jj]*moments[qinds[ii][jj]];
				sxv[i][j] /= moments[0];
			}
		}
		if(grid->fields[PSI_GRID_VV]) {
			for(i = 0; i < PSI_NDIM; ++i) 
			for(j = 0; j < PSI_NDIM; ++j)  {
				svv[i][j] = moments[0]*(vel[0].xyz[i]-cv[i])*(vel[0].xyz[j]-cv[j]);
				for(ii = 0; ii < PSI_NDIM; ++ii) {
					svv[i][j] += DV[i][ii]*moments[1+ii]*(vel[0].xyz[j]-cv[j]);
					svv[i][j] += DV[j][ii]*moments[1+ii]*(vel[0].xyz[i]-cv[i]);
				}
				for(ii = 0; ii < PSI_NDIM; ++ii)
				for(jj = 0; jj < PSI_NDIM; ++jj)
					svv[i][j] += DV[i][ii]*DV[j][jj]*moments[qinds[ii][jj]];
				svv[i][j] /= moments[0];
			}
		}

		// update everything online 
		// TODO: go back to adding up just the raw moments...
		cwght = mass*moments[0]/(grid->fields[PSI_GRID_M][gridind] + mass*moments[0]); 
		dwght = grid->fields[PSI_GRID_M][gridind]/(grid->fields[PSI_GRID_M][gridind] + mass*moments[0]); 
		if(grid->fields[PSI_GRID_XX]) for(i = 0; i < PSI_NDIM; ++i) for(j = 0; j < PSI_NDIM; ++j) {
			grid->fields[PSI_GRID_XX][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] = 
				dwght*grid->fields[PSI_GRID_XX][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] + cwght*sxx[i][j] 
				+ cwght*dwght*(cm[i]-grid->fields[PSI_GRID_X][PSI_NDIM*gridind+i])*(cm[j]-grid->fields[PSI_GRID_X][PSI_NDIM*gridind+j]);
		}
		if(grid->fields[PSI_GRID_XV]) for(i = 0; i < PSI_NDIM; ++i) for(j = 0; j < PSI_NDIM; ++j) {
			grid->fields[PSI_GRID_XV][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] = 
				dwght*grid->fields[PSI_GRID_XV][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] + cwght*sxv[i][j] 
				+ cwght*dwght*(cm[i]-grid->fields[PSI_GRID_X][PSI_NDIM*gridind+i])*(cv[j]-grid->fields[PSI_GRID_V][PSI_NDIM*gridind+j]);
		}
		if(grid->fields[PSI_GRID_VV]) for(i = 0; i < PSI_NDIM; ++i) for(j = 0; j < PSI_NDIM; ++j) {
			grid->fields[PSI_GRID_VV][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] = 
				dwght*grid->fields[PSI_GRID_VV][PSI_NDIM*PSI_NDIM*gridind+PSI_NDIM*i+j] + cwght*svv[i][j] 
				+ cwght*dwght*(cv[i]-grid->fields[PSI_GRID_V][PSI_NDIM*gridind+i])*(cv[j]-grid->fields[PSI_GRID_V][PSI_NDIM*gridind+j]);
		}
		if(grid->fields[PSI_GRID_X]) for(i = 0; i < PSI_NDIM; ++i) 
			grid->fields[PSI_GRID_X][PSI_NDIM*gridind+i] += cwght*(cm[i]-grid->fields[PSI_GRID_X][PSI_NDIM*gridind+i]);
		if(grid->fields[PSI_GRID_V]) for(i = 0; i < PSI_NDIM; ++i) 
			grid->fields[PSI_GRID_V][PSI_NDIM*gridind+i] += cwght*(cv[i]-grid->fields[PSI_GRID_V][PSI_NDIM*gridind+i]);

		if(grid->fields[PSI_GRID_M]) grid->fields[PSI_GRID_M][gridind] += mass*moments[0];
	}
}

****/
