template<typename part, typename part_info, typename part_dataType>
void test_id(Particles<part, part_info, part_dataType> * pcls, long *ID, double boxsize)
{
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
	psi_rvec gpos[(1<<meshdim)*(meshdim+1)];
	psi_rvec grbox[(1<<meshdim)*2];

	// constant-size buffer for max refinement level
	psi_rtree_query qry;
	psi_tet_buffer tetbuf, tetbuf1;
	psi_tet_buffer_init(&tetbuf, 0, 0);

    double window[2][3];
    window[0][0] = 0.;
    window[0][1] = 0.;
    window[0][2] = 0.;
    window[1][0] = 256.;
    window[1][1] = 256.;
    window[1][2] = 256.;

    Site xPart(pcls->lattice());
    typename std::list<part>::iterator it;

    double pos[8][3], vel[8][3];

    for (xPart.first(); xPart.test(); xPart.next()){
        for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it){
            for (psi_int i = 0; i < 8; ++i){
                // find the right 8 particles
                if (it->ID == ID[i]){
                    pos[i][0] = it->pos[0];
                    pos[i][1] = it->pos[1];
                    pos[i][2] = it->pos[2];
                    vel[i][0] = it->vel[0];
                    vel[i][1] = it->vel[1];
                    vel[i][2] = it->vel[2];
                }
            }
        }
    }
    for (xPart.first(); xPart.test(); xPart.next()){
        for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it){

            // a local copy of the element, in case it's modified
            // make it periodic, get its bounding box, and check it against the grid
            tmass = 1;//mesh->mass[e]; 
            if(!psi_aabb_periodic(pos, trbox, grid->window, mesh)){
                continue;
            }

            // refine elements into the tet buffer 
            psi_tet_buffer_refine(&tetbuf, tpos, tvel, tmass, 8);

            // linear in the density
            if(mode == PSI_MODE_DENSITY){

                // loop over each tet in the buffer
                // copy the position to the ghost array and compute its aabb
                // make ghosts and sample tets to the grid
                for(psi_int t = 0; t < tetbuf.num; ++t){
                    memcpy(gpos, &tetbuf.pos[(meshdim+1)*t], (meshdim+1)*sizeof(psi_rvec));
                    psi_aabb(gpos, meshdim+1,grbox);
                    //psi_make_ghosts(gpos, grbox, &nghosts, (meshdim+1), grid->window, mesh);
                    nghosts = 1;
                    for(psi_int g = 0; g < nghosts; ++g){
                        psi_voxelize_tet(&gpos[(meshdim+1)*g], &tetbuf.vel[(meshdim+1)*t], tetbuf.mass[t], &grbox[2*g], grid);
                    }
                }
            }
        }
    }
}

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


