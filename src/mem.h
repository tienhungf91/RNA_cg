#ifndef MEM_H
#define MEM_H

//////////////////////////////////////////////////////////////////
void allocate_mem (_topol_struct           &topol_struct,
                   _crowder_struct         &crowder_struct,
                   _coord_vel_force_struct &coord_vel_force_struct) {
/////// Remember that crowder indices start at 0, instead of 1 as in RNA

//    crowder_struct.crowder_key = (int *) calloc (topol_struct.Nres + 1, sizeof(int));
    crowder_struct.crowder_key.resize (topol_struct.Nres + 1, 0);

    crowder_struct.crowder_mass.resize    (N_CRWD, 0);
    crowder_struct.crowder_content.resize (N_CRWD);

//    for (int i = 0; i < N_CRWD; i++) {
    for (int i = 0; i < 1; i++) {
        crowder_struct.crowder_mass [i] = crowder_struct.N_crwd [i];
        crowder_struct.crowder_content [i].push_back (0);     // Not using "0" index for crowder_content
//        crowder_struct.crowder_content [i] = (int *) calloc (crowder_struct.crowder_mass [i] + 1, sizeof(int));
    }

    crowder_struct.crowder_key.resize (topol_struct.Nres + 1, -1);
    topol_struct.Nsite.resize         (topol_struct.Nres + 1);

//    crowder_struct.Nsite = (int *) calloc (crowder_struct.Nsite + 1, sizeof(int));
//    topol_struct.Nsite = (int *) realloc (topol_struct.Nsite, (topol_struct.Nres + 1) * sizeof(int));
/*    topol_struct.RIGID_SET.resize     (topol_struct.Nres + 1);
    topol_struct.RIGID_END.resize     (topol_struct.Nres + 1);
    topol_struct.RX.resize            (topol_struct.Nres + 1);
    topol_struct.RY.resize            (topol_struct.Nres + 1);
    topol_struct.RZ.resize            (topol_struct.Nres + 1);*/

    topol_struct.atom_key = (int **) realloc (topol_struct.atom_key, (topol_struct.Nres + 1) * sizeof(int *));
    topol_struct.Natm = topol_struct.Natm_biomol;    // start here and update later

    for (int i = 1; i <= topol_struct.Nres - topol_struct.Nres_biomol; i++) {
        int index = i + topol_struct.Nres_biomol;
        int type;
        if      (i <= crowder_struct.N_crwd [0])                             type = 0; /////// MG
        else if (i <= crowder_struct.N_crwd [0] + crowder_struct.N_crwd [1]) type = 1; /////// Cl
        else                                                                 type = 2; /////// K

        crowder_struct.crowder_content [type].push_back (index);
        crowder_struct.crowder_key [index] = type;
        topol_struct.Nsite         [index] = crowder_struct.Nsite [type];
//        topol_struct.RIGID_SET     [index] = crowder_struct.RIGID_SET_crwd [crowder_type];
//        topol_struct.RIGID_END     [index] = crowder_struct.RIGID_END_crwd [crowder_type];

//        topol_struct.atom_key [index] = NULL;
        topol_struct.atom_key [index] = (int *) calloc (topol_struct.Nsite [index], sizeof(int));
        for (int j = 0; j < topol_struct.Nsite [index]; j++) {
            topol_struct.Natm ++;
            topol_struct.atom_key [index][j] = topol_struct.Natm;

            topol_struct.part_key.push_back (crowder_struct.part_key_crwd [type][j]);
            topol_struct.maxi_key.push_back (crowder_struct.maxi_key_crwd [type][j]);

/*            topol_struct.part_key = (int *) realloc (topol_struct.part_key, (topol_struct.Natm + 1) * sizeof(int));
            topol_struct.part_key [topol_struct.Natm] = crowder_struct.part_key_crwd [crowder_type][j];

            topol_struct.maxi_key = (int *) realloc (topol_struct.maxi_key, (topol_struct.Natm + 1) * sizeof(int));
            topol_struct.maxi_key [topol_struct.Natm] = crowder_struct.maxi_key_crwd [crowder_type][j];

            topol_struct.INDX = (int *) realloc (topol_struct.INDX, (topol_struct.Natm + 1) * sizeof(int));
            topol_struct.JNDX = (int *) realloc (topol_struct.JNDX, (topol_struct.Natm + 1) * sizeof(int));*/

//            topol_struct.INDX.push_back (index);
//            topol_struct.JNDX.push_back (j);
//            topol_struct.INDX [topol_struct.Natm] = i;
//            topol_struct.JNDX [topol_struct.Natm] = j;

/*            topol_struct.RX [index].resize (topol_struct.Nsite [index] + 1, 0);
            topol_struct.RY [index].resize (topol_struct.Nsite [index] + 1, 0);
            topol_struct.RZ [index].resize (topol_struct.Nsite [index] + 1, 0);*/

            //coord_vel_force_struct.coordx = (double *) realloc (coord_vel_force_struct.coordx, (topol_struct.Natm + 1) * sizeof(double));
            //coord_vel_force_struct.coordy = (double *) realloc (coord_vel_force_struct.coordy, (topol_struct.Natm + 1) * sizeof(double));
            //coord_vel_force_struct.coordz = (double *) realloc (coord_vel_force_struct.coordz, (topol_struct.Natm + 1) * sizeof(double));

            //coord_vel_force_struct.coordx [topol_struct.Natm] = 0;
            //coord_vel_force_struct.coordy [topol_struct.Natm] = 0;
            //coord_vel_force_struct.coordz [topol_struct.Natm] = 0;
        }
    }

    coord_vel_force_struct.coordx.resize (topol_struct.Natm + 1, 0);
    coord_vel_force_struct.coordy.resize (topol_struct.Natm + 1, 0);
    coord_vel_force_struct.coordz.resize (topol_struct.Natm + 1, 0);

    coord_vel_force_struct.coord_oldx.resize (topol_struct.Natm + 1, 0);
    coord_vel_force_struct.coord_oldy.resize (topol_struct.Natm + 1, 0);
    coord_vel_force_struct.coord_oldz.resize (topol_struct.Natm + 1, 0);

    coord_vel_force_struct.velox.resize (topol_struct.Natm + 1, 0);
    coord_vel_force_struct.veloy.resize (topol_struct.Natm + 1, 0);
    coord_vel_force_struct.veloz.resize (topol_struct.Natm + 1, 0);

/*    coord_vel_force_struct.coordx0 = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.coordy0 = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.coordz0 = (double *) calloc (topol_struct.Natm + 1, sizeof(double));

    coord_vel_force_struct.coord_oldx = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.coord_oldy = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.coord_oldz = (double *) calloc (topol_struct.Natm + 1, sizeof(double));

    coord_vel_force_struct.velox      = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.veloy      = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.veloz      = (double *) calloc (topol_struct.Natm + 1, sizeof(double));*/

    coord_vel_force_struct.forcex     = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.forcey     = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.forcez     = (double *) calloc (topol_struct.Natm + 1, sizeof(double));

    coord_vel_force_struct.flx        = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.fly        = (double *) calloc (topol_struct.Natm + 1, sizeof(double));
    coord_vel_force_struct.flz        = (double *) calloc (topol_struct.Natm + 1, sizeof(double));

    // Make it a multiple of 6
    int tmp = topol_struct.Natm % 6;
    coord_vel_force_struct.maxwell_force.resize (topol_struct.Natm - tmp + 7, 0);
}

////////////////////////////////////////////////////////////////
void free_hydrogen_bonds (_topol_struct   &topol_struct,
                          _hydbond_struct &hydbond_struct) {
//    free (hydbond_struct.hb_dode);
    hydbond_struct.hb_dode.clear ();

    for (int i = 1; i <= hydbond_struct.hb_N; i++) free (hydbond_struct.hb_code [i]);
    free (hydbond_struct.hb_code);

//    free (hydbond_struct.hb_RES1); 
//    free (hydbond_struct.hb_RES2);

//    free (hydbond_struct.ATOM_HB_N);
//    hydbond_struct.hb_RES1.clear ();
//    hydbond_struct.hb_RES2.clear ();
    hydbond_struct.ATOM_HB_N.clear ();

    for (int i = 1; i <= hydbond_struct.hb_N; i++) free (hydbond_struct.ATOM_HB [i]);
    free (hydbond_struct.ATOM_HB);

/*    free (hydbond_struct.HB_EXCESS);
    free (hydbond_struct.EXCESS);
    free (hydbond_struct.HB_A);

    free (hydbond_struct.HB_ATOM_N);*/

    hydbond_struct.HB_EXCESS.clear ();
    hydbond_struct.EXCESS.clear ();
    hydbond_struct.HB_A.clear ();
    hydbond_struct.HB_ATOM_N.clear ();

    for (int i = 1; i <= hydbond_struct.HB_Natm; i++) free (hydbond_struct.HB_ATOM [i]);
    free (hydbond_struct.HB_ATOM);

//    free (hydbond_struct.HB_PAIR_N);
    hydbond_struct.HB_PAIR_N.clear ();

    for (int i = 1; i <= topol_struct.Npair_biomol; i++) free (hydbond_struct.HB_PAIR [i]);
    free (hydbond_struct.HB_PAIR);

//    hydbond_struct.HB_Natm = 0;
//    free (hydbond_struct.HB_Nsite); hydbond_struct.HB_Nsite = NULL;
    hydbond_struct.HB_Nsite.clear ();

    for (int i = 1; i <= topol_struct.Nres_biomol; i++) free (hydbond_struct.HB_atom_key [i]);
    free (hydbond_struct.HB_atom_key); hydbond_struct.HB_atom_key = NULL;
//    hydbond_struct.HB_atom_key.clear ();

/*    free (hydbond_struct.HB_INDX); hydbond_struct.HB_INDX = NULL; 
    free (hydbond_struct.HB_JNDX); hydbond_struct.HB_JNDX = NULL;

    free (hydbond_struct.VALENCE); hydbond_struct.VALENCE = NULL;

    hydbond_struct.hb_N = 0;

    free (hydbond_struct.hb_k10);    hydbond_struct.hb_k10 = NULL;
    free (hydbond_struct.hb_k11);    hydbond_struct.hb_k11 = NULL;
    free (hydbond_struct.hb_k12);    hydbond_struct.hb_k12 = NULL;
    free (hydbond_struct.hb_k20);    hydbond_struct.hb_k20 = NULL;
    free (hydbond_struct.hb_k21);    hydbond_struct.hb_k21 = NULL;
    free (hydbond_struct.hb_k22);    hydbond_struct.hb_k22 = NULL;
    free (hydbond_struct.hb_r);      hydbond_struct.hb_r = NULL;
    free (hydbond_struct.hb_theta1); hydbond_struct.hb_theta1 = NULL;
    free (hydbond_struct.hb_theta2); hydbond_struct.hb_theta2 = NULL;
    free (hydbond_struct.hb_psi);    hydbond_struct.hb_psi = NULL;
    free (hydbond_struct.hb_psi1);   hydbond_struct.hb_psi1 = NULL;
    free (hydbond_struct.hb_psi2);   hydbond_struct.hb_psi2 = NULL;
    free (hydbond_struct.hb_E);      hydbond_struct.hb_E = NULL;
    free (hydbond_struct.hb_psi0);   hydbond_struct.hb_psi0 = NULL;
    free (hydbond_struct.hb_psi10);  hydbond_struct.hb_psi10 = NULL;
    free (hydbond_struct.hb_psi20);  hydbond_struct.hb_psi20 = NULL;
    free (hydbond_struct.hb_energy); hydbond_struct.hb_energy = NULL;
    free (hydbond_struct.hb_status); hydbond_struct.hb_status = NULL;
    free (hydbond_struct.hb_K);      hydbond_struct.hb_K = NULL;*/

//    hydbond_struct.HB_INDX.clear ();    hydbond_struct.HB_JNDX.clear ();
    hydbond_struct.VALENCE.clear ();
    hydbond_struct.hb_k10.clear ();     hydbond_struct.hb_k11.clear ();     hydbond_struct.hb_k12.clear ();
    hydbond_struct.hb_k20.clear ();     hydbond_struct.hb_k21.clear ();     hydbond_struct.hb_k22.clear ();
    hydbond_struct.hb_r.clear ();       hydbond_struct.hb_theta1.clear ();  hydbond_struct.hb_theta2.clear ();
    hydbond_struct.hb_psi.clear ();     hydbond_struct.hb_psi1.clear ();    hydbond_struct.hb_psi2.clear ();
    hydbond_struct.hb_psi0.clear ();    hydbond_struct.hb_psi10.clear ();   hydbond_struct.hb_psi20.clear ();
    hydbond_struct.hb_E.clear ();       hydbond_struct.hb_energy.clear ();  hydbond_struct.hb_status.clear ();
    hydbond_struct.hb_K.clear ();
}

/////////////////////////////////////////////////////////////////////////
void free_stacks (//_topol_struct &topol_struct,
                  _stack_struct &stack_struct) {
//    for (int k = 1; k < stack_struct.st_N + 1; k++) free (stack_struct.ATOM_ST [k]);
//    free (stack_struct.ATOM_ST);
//    stack_struct.ATOM_ST.clear();

//    free (stack_struct.ST_EXCESS);
//    stack_struct.ST_EXCESS.clear ();

//    for (int k = 1; k <= topol_struct.Nres_biomol; k++) free (stack_struct.ST_ATOM [k]);
//    free (stack_struct.ST_ATOM);
//    free (stack_struct.ST_ATOM_N);
    stack_struct.ST_ATOM_N.clear ();

/*    stack_struct.s3_N = 0; stack_struct.st_N = 0;

    free (stack_struct.st_i);      stack_struct.st_i = NULL;
    free (stack_struct.st_j);      stack_struct.st_j = NULL;
    free (stack_struct.st_r);      stack_struct.st_r = NULL;

    free (stack_struct.st_theta1); stack_struct.st_theta1 = NULL;
    free (stack_struct.st_theta2); stack_struct.st_theta2 = NULL;
    free (stack_struct.st_psi);    stack_struct.st_psi = NULL;
    free (stack_struct.st_psi1);   stack_struct.st_psi1 = NULL;
    free (stack_struct.st_psi2);   stack_struct.st_psi2 = NULL;

    free (stack_struct.st_E);      stack_struct.st_E = NULL;

    free (stack_struct.st_psi0);   stack_struct.st_psi0 = NULL;
    free (stack_struct.st_psi10);  stack_struct.st_psi10 = NULL;
    free (stack_struct.st_psi20);  stack_struct.st_psi20 = NULL;

    free (stack_struct.st_r2);     stack_struct.st_r2 = NULL;

    free (stack_struct.st_energy); stack_struct.st_energy = NULL;
    free (stack_struct.st_status); stack_struct.st_status = NULL;
    free (stack_struct.e3_stack);  stack_struct.e3_stack = NULL;*/

    stack_struct.st_i2.clear ();   stack_struct.st_j2.clear ();      stack_struct.st_r.clear ();
    stack_struct.st_i3.clear ();   stack_struct.st_j3.clear ();
    stack_struct.st_E.clear ();    stack_struct.st_theta1.clear ();  stack_struct.st_theta2.clear ();
    stack_struct.st_psi.clear ();  stack_struct.st_psi1.clear ();    stack_struct.st_psi2.clear ();
    stack_struct.st_psi0.clear (); stack_struct.st_psi10.clear ();   stack_struct.st_psi20.clear ();
    stack_struct.st_energy.clear ();  stack_struct.st_status.clear ();
    stack_struct.e3_stack.clear ();
}

#endif
