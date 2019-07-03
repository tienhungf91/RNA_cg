#ifndef CROWDER_H
#define CROWDER_H

void list_crowder_types (int ion_type,
                         _crowder_struct &crowder_struct) {

//    crowder_struct.crowder_mass    = (int *)  calloc (N_CRWD, sizeof(int));
//    crowder_struct.crowder_content = (int **) calloc (N_CRWD, sizeof(int *));
    crowder_struct.part_key_crwd.resize (N_CRWD);
    crowder_struct.maxi_key_crwd.resize (N_CRWD);

    for (int i = 0; i < N_CRWD; i++)
        switch (i) {
            case 0:   //Mg
                crowder_struct.Nsite [i] = 1;
//                crowder_struct.RIGID_SET_crwd [i] = 0;
//                crowder_struct.RIGID_END_crwd [i] = -1;

//                crowder_struct.part_key_crwd [i] = NULL;
//                crowder_struct.part_key_crwd [i] = (int *) calloc (crowder_struct.Nsite [i], sizeof(int));
//                crowder_struct.part_key_crwd [i][0] = 3;
                crowder_struct.part_key_crwd [i].push_back (3);

//                crowder_struct.maxi_key_crwd [i] = NULL;
//                crowder_struct.maxi_key_crwd [i] = (int *) calloc (crowder_struct.Nsite [i], sizeof(int));
//                crowder_struct.maxi_key_crwd [i][0] = 7;
                crowder_struct.maxi_key_crwd [i].push_back (6 + ion_type);

/*                crowder_struct.RX_crwd [i] = NULL;
                crowder_struct.RY_crwd [i] = NULL;
                crowder_struct.RZ_crwd [i] = NULL;

                crowder_struct.RX_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));
                crowder_struct.RY_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));
                crowder_struct.RZ_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));

                crowder_struct.RX_crwd [i][0] = 0;
                crowder_struct.RY_crwd [i][0] = 0;
                crowder_struct.RZ_crwd [i][0] = 0;*/
                break;

            case 1://Cl
                crowder_struct.Nsite [i] = 1;
//                crowder_struct.RIGID_SET_crwd [i] = 0;
//                crowder_struct.RIGID_END_crwd [i] = -1;

//                crowder_struct.part_key_crwd [i] = NULL;
//                crowder_struct.part_key_crwd [i] = (int *) calloc (crowder_struct.Nsite [i], sizeof(int));
//                crowder_struct.part_key_crwd [i][0] = 4;
                crowder_struct.part_key_crwd [i].push_back (4);

//                crowder_struct.maxi_key_crwd [i] = NULL;
//                crowder_struct.maxi_key_crwd [i] = (int *) calloc (crowder_struct.Nsite [i], sizeof(int));
//                crowder_struct.maxi_key_crwd [i][0] = 8;
                crowder_struct.maxi_key_crwd [i].push_back (8);

/*                crowder_struct.RX_crwd [i] = NULL;
                crowder_struct.RY_crwd [i] = NULL;
                crowder_struct.RZ_crwd [i] = NULL;

                crowder_struct.RX_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));
                crowder_struct.RY_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));
                crowder_struct.RZ_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));

                crowder_struct.RX_crwd [i][0] = 0;
                crowder_struct.RY_crwd [i][0] = 0;
                crowder_struct.RZ_crwd [i][0] = 0;*/
                break;

            case 2://K
                crowder_struct.Nsite [i] = 1;
//                crowder_struct.RIGID_SET_crwd [i] = 0;
//                crowder_struct.RIGID_END_crwd [i] = -1;

//                crowder_struct.part_key_crwd [i] = NULL;
//                crowder_struct.part_key_crwd [i] = (int *) calloc (crowder_struct.Nsite [i], sizeof(int));
//                crowder_struct.part_key_crwd [i][0] = 5;
                crowder_struct.part_key_crwd [i].push_back (5);

//                crowder_struct.maxi_key_crwd [i] = NULL;
//                crowder_struct.maxi_key_crwd [i] = (int *) calloc (crowder_struct.Nsite [i], sizeof(int));
//                crowder_struct.maxi_key_crwd [i][0] = 9;
                crowder_struct.maxi_key_crwd [i].push_back (9);

/*                crowder_struct.RX_crwd [i] = NULL;
                crowder_struct.RY_crwd [i] = NULL;
                crowder_struct.RZ_crwd [i] = NULL;

                crowder_struct.RX_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));
                crowder_struct.RY_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));
                crowder_struct.RZ_crwd [i] = (double *) calloc (crowder_struct.Nsite [i], sizeof(double));

                crowder_struct.RX_crwd [i][0] = 0;
                crowder_struct.RY_crwd [i][0] = 0;
                crowder_struct.RZ_crwd [i][0] = 0;*/
                break;

            default:
                break;
        }
}

/////////////////////////////////////////////////////////////////////////////////////////////
void add_crowder (//int tp,              //// tp: type of crowders
                  int index, double * box,
                  const _topol_struct     &topol_struct,
                  _coord_vel_force_struct &coord_vel_force_struct,
                  bool ion_close, double mean_dist) {
    bool added = 0, added_close = 0;
    while (not added && not added_close) {
        added = 1;

//#ifdef _OPENMP
//        unsigned int myseed = omp_get_thread_num();
//
//        double a1 = (double) rand_r (&myseed) / RAND_MAX;
//        double tmp1 = a1 * box[0];
//
//        a1 = (double) rand_r (&myseed) / RAND_MAX;
//        double tmp2 = a1 * box[1];
//
//        a1 = (double) rand_r (&myseed) / RAND_MAX;
//        double tmp3 = a1 * box[2];
//#else
        double a1 = (double) rand () / RAND_MAX;
        double tmp1 = a1 * box[0];

        a1 = (double) rand () / RAND_MAX;
        double tmp2 = a1 * box[1];

        a1 = (double) rand () / RAND_MAX;
        double tmp3 = a1 * box[2];
//#endif

/*        _Dcoordinate n1;
        n1.x = rand ();
        a1 = (double) rand () / RAND_MAX;
        if (a1 > 0.5) n1.x = -n1.x; //n1 [0] = - n1 [0];

        n1.y = rand ();
        a1 = (double) rand () / RAND_MAX;
        if (a1 > 0.5) n1.y = -n1.y; //n1 [1] = - n1 [1];

        n1.z = rand ();
        a1 = (double) rand () / RAND_MAX;
        if (a1 > 0.5) n1.z = -n1.z; //n1 [2] = - n1 [2];

        a1 = sqrt (n1.x*n1.x + n1.y*n1.y + n1.z*n1.z);

        n1.x /= a1;       n1.y /= a1;       n1.z /= a1;

        _Dcoordinate n2 = gen_perpendicular_vector (n1);
        _Dcoordinate n3 = cross_product (n1, n2);*/

        for (int j = 0; j < topol_struct.Nsite [index]; j++) {
            if (not added || added_close) break;
            int k = topol_struct.atom_key [index][j];
            int In1 = topol_struct.maxi_key [k];

            coord_vel_force_struct.coordx [k] = tmp1;
            coord_vel_force_struct.coordy [k] = tmp2;
            coord_vel_force_struct.coordz [k] = tmp3;

            ////// check for clashing with RNA and other ions
            //        for (int i2 = 1; i2 <= topol_struct.Nres_biomol; i2++)
            for (int i2 = 1; i2 < index; i2++) {
                if (not added) break;
//                for (int j2 = 0; j2 < topol_struct.Nsite [i2]; j2++) {
                for (int j2 = 0; j2 < 3; j2++) {
                    if (not added) break;
                    int k2 = topol_struct.atom_key [i2][j2];
                    if (k2 == 0) continue;
                    int In2 = topol_struct.maxi_key [k2];

                    int In12;
                    if (In1 < In2)
                        In12 = (In1 - 1) * NA - (In1 - 1) * In1 / 2 + In2;

                    else In12 = (In2 - 1) * NA - (In2 - 1) * In2 / 2 + In1;

                    double dist [3];
                    dist[0] = coord_vel_force_struct.coordx [k] - coord_vel_force_struct.coordx [k2];
                    dist[1] = coord_vel_force_struct.coordy [k] - coord_vel_force_struct.coordy [k2];
                    dist[2] = coord_vel_force_struct.coordz [k] - coord_vel_force_struct.coordz [k2];

                    half_shift (dist, box);

                    double cutoff_dist_sqr = 0;
                    double mean_d2 = mean_dist * mean_dist;

                    // Add ions normally
                    if (not ion_close) {
                        // if (i2 <= topol_struct.Nres_biomol)
                        //     cutoff_dist_sqr = std::min (100., mean_d2); // Keep ions away from RNA
                        // else if (In2 != In1)
                        if (In2 != In1)
                            //if (In1 == 8 || In2 == 8) cutoff_dist_sqr = std::min (64., 0.5*mean_d2); // Keep Cl- away from others
                            //else
                            cutoff_dist_sqr = topol_struct.D2_LJ_OVERLAP [In12];
                        else cutoff_dist_sqr = mean_d2;

                        if (dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2] < cutoff_dist_sqr)
                            added = 0;
                    } else {
                        // Keep ions close to the biomolecule
                        double d2 = dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
                        if (d2 < mean_d2 && d2 > 100.)
                            added_close = 1;
                    }
                }
            }
        }
    }
}

#endif
