#ifndef MANI_COORD_H
#define MANI_COORD_H

////////////////////////////////////////////////////////////////////////////
void set_box (double * box,
              const _topol_struct     &topol_struct,
              _coord_vel_force_struct &coord_vel_force_struct) {

    ///////////// RNA axes /////////////
    double a0 = 0, a1 = 0, a2 = 0;
    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        a0 += coord_vel_force_struct.coordx [i];
        a1 += coord_vel_force_struct.coordy [i];
        a2 += coord_vel_force_struct.coordz [i];
    }

    a0 /= topol_struct.Natm_biomol;
    a1 /= topol_struct.Natm_biomol;
    a2 /= topol_struct.Natm_biomol;

    double A = 0, B = 0, C = 0, F = 0, G = 0, H = 0; 

    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        double x = coord_vel_force_struct.coordx [i] - a0;
        double y = coord_vel_force_struct.coordy [i] - a1;
        double z = coord_vel_force_struct.coordz [i] - a2;

        A += y*y + z*z;
        B += x*x + z*z;
        C += x*x + y*y;
        F += y*z;
        G += x*z;
        H += x*y;
    }

    _Dcoordinate ep1, ep2, ep3;
    compute_ep_vector (A, B, C, F, G, H, ep1, ep2, ep3);

    /// Rotation
    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        double x_new = coord_vel_force_struct.coordx [i] * ep1.x + coord_vel_force_struct.coordy [i] * ep1.y + coord_vel_force_struct.coordz [i] * ep1.z;
        double y_new = coord_vel_force_struct.coordx [i] * ep2.x + coord_vel_force_struct.coordy [i] * ep2.y + coord_vel_force_struct.coordz [i] * ep2.z;
        double z_new = coord_vel_force_struct.coordx [i] * ep3.x + coord_vel_force_struct.coordy [i] * ep3.y + coord_vel_force_struct.coordz [i] * ep3.z;

        coord_vel_force_struct.coordx [i] = x_new;
        coord_vel_force_struct.coordy [i] = y_new;
        coord_vel_force_struct.coordz [i] = z_new;
    }

    ////////// Move the RNA to the box center
    double x = 0.5*box[0] - coord_vel_force_struct.coordx [topol_struct.bead_center];
    double y = 0.5*box[1] - coord_vel_force_struct.coordy [topol_struct.bead_center];
    double z = 0.5*box[2] - coord_vel_force_struct.coordz [topol_struct.bead_center];

    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        coord_vel_force_struct.coordx [i] += x;
        coord_vel_force_struct.coordy [i] += y;
        coord_vel_force_struct.coordz [i] += z;
    }
}

///////////////////////////////////////////////////////////////////////////
void generate_atom_velocity (int i, int j, double T,     // temperature
                             const _topol_struct     &topol_struct,
                             _coord_vel_force_struct &coord_vel_force_struct) {
    double a1, a2, a3, a4, a5, a6;

    int k = topol_struct.atom_key [i][j];
    if (k == 0) return;
    int l = topol_struct.maxi_key [k];

//#ifdef _OPENMP
//    unsigned int myseed = omp_get_thread_num() * T * 1000;
//
//    do a1 = (double) rand_r (&myseed) / RAND_MAX; while (a1 < 1e-6);
//    a2 = (double) rand_r (&myseed) / RAND_MAX;
//
//    do a3 = (double) rand_r (&myseed) / RAND_MAX; while (a3 < 1e-6);
//    a4 = (double) rand_r (&myseed) / RAND_MAX;
//
//    do a5 = (double) rand_r (&myseed) / RAND_MAX; while (a5 < 1e-6);
//    a6 = (double) rand_r (&myseed) / RAND_MAX;
//#else
    do a1 = (double) rand () / RAND_MAX; while (a1 < 1e-6);
    a2 = (double) rand () / RAND_MAX;

    do a3 = (double) rand () / RAND_MAX; while (a3 < 1e-6);
    a4 = (double) rand () / RAND_MAX;

    do a5 = (double) rand () / RAND_MAX; while (a5 < 1e-6);
    a6 = (double) rand () / RAND_MAX;
//#endif

    double tmp1 = sqrt (T / topol_struct.MASS [l]) * sqrt (-2*log (a1));
    double tmp2 = 2*PI * a2;
    coord_vel_force_struct.velox [k] = tmp1 * cos (tmp2);

    tmp1 = sqrt (T / topol_struct.MASS [l]) * sqrt (-2*log (a3));
    tmp2 = 2*PI * a4;
    coord_vel_force_struct.veloy [k] = tmp1 * cos (tmp2);

    tmp1 = sqrt (T / topol_struct.MASS [l]) * sqrt (-2*log (a5));
    tmp2 = 2*PI * a6;
    coord_vel_force_struct.veloz [k] = tmp1 * cos (tmp2);
}

////////////////////////////////////////////////////////////////////////////
void gen_old_coord (int Natm, double dt,
                    _coord_vel_force_struct &coord_vel_force_struct) {
    for (int i = 1; i <= Natm; i++) {
        coord_vel_force_struct.coord_oldx [i] = coord_vel_force_struct.coordx [i] - coord_vel_force_struct.velox [i] * dt;
        coord_vel_force_struct.coord_oldy [i] = coord_vel_force_struct.coordy [i] - coord_vel_force_struct.veloy [i] * dt;
        coord_vel_force_struct.coord_oldz [i] = coord_vel_force_struct.coordz [i] - coord_vel_force_struct.veloz [i] * dt;
    }
}

////////////////////////////////////////////////////////////////
void leap_frog_move_atom (int k, double temp, double dt,
//                          long int step, FILE * f_error,
                          const _topol_struct     &topol_struct,
                          _coord_vel_force_struct &coord_vel_force_struct) {
    int l = topol_struct.maxi_key [k];

    double x_tmp  = coord_vel_force_struct.coordx [k];
    double y_tmp  = coord_vel_force_struct.coordy [k];
    double z_tmp  = coord_vel_force_struct.coordz [k];
    double velox = coord_vel_force_struct.velox [k];
    double veloy = coord_vel_force_struct.veloy [k];
    double veloz = coord_vel_force_struct.veloz [k];

    coord_vel_force_struct.velox [k] = topol_struct.K1 [l] * velox + topol_struct.K2 [l] * coord_vel_force_struct.forcex [k];
    coord_vel_force_struct.veloy [k] = topol_struct.K1 [l] * veloy + topol_struct.K2 [l] * coord_vel_force_struct.forcey [k];
    coord_vel_force_struct.veloz [k] = topol_struct.K1 [l] * veloz + topol_struct.K2 [l] * coord_vel_force_struct.forcez [k];

    /// velocity scaling
    double vsq = topol_struct.MASS [l] * (coord_vel_force_struct.velox [k] * coord_vel_force_struct.velox [k] +\
                                          coord_vel_force_struct.veloy [k] * coord_vel_force_struct.veloy [k] +\
                                          coord_vel_force_struct.veloz [k] * coord_vel_force_struct.veloz [k]) / (300.0*temp);
    if (vsq > 1.0) {
//        fprintf (f_error, "step = %10ld         atom_key = %6d         vsq = %8.3e\n", step, k, vsq);
        vsq = 1.0 / sqrt (vsq);
        coord_vel_force_struct.velox [k] *= vsq;
        coord_vel_force_struct.veloy [k] *= vsq;
        coord_vel_force_struct.veloz [k] *= vsq;
    }

    coord_vel_force_struct.coordx [k] = 2*x_tmp - coord_vel_force_struct.coord_oldx [k] + (coord_vel_force_struct.velox [k] - velox) * dt;
    coord_vel_force_struct.coordy [k] = 2*y_tmp - coord_vel_force_struct.coord_oldy [k] + (coord_vel_force_struct.veloy [k] - veloy) * dt;
    coord_vel_force_struct.coordz [k] = 2*z_tmp - coord_vel_force_struct.coord_oldz [k] + (coord_vel_force_struct.veloz [k] - veloz) * dt;

    coord_vel_force_struct.coord_oldx [k] = x_tmp;
    coord_vel_force_struct.coord_oldy [k] = y_tmp;
    coord_vel_force_struct.coord_oldz [k] = z_tmp;
}

//////////////////////////////////////////////////////////////////
void adjust_box (int bead_center, int Natm,
                 double * box,
                 _coord_vel_force_struct &coord_vel_force_struct) {

    double dx = 0.5*box[0] - coord_vel_force_struct.coordx [bead_center];
    double dy = 0.5*box[1] - coord_vel_force_struct.coordy [bead_center];
    double dz = 0.5*box[2] - coord_vel_force_struct.coordz [bead_center];

    for (int k = 1; k <= Natm; k++) {
        coord_vel_force_struct.coordx [k]     += dx;
        coord_vel_force_struct.coord_oldx [k] += dx;

        coord_vel_force_struct.coordy [k]     += dy;
        coord_vel_force_struct.coord_oldy [k] += dy;

        coord_vel_force_struct.coordz [k]     += dz;
        coord_vel_force_struct.coord_oldz [k] += dz;
    }
}

//////////////////////////////////////////////////////////////////////////////////
void wrap_coord (const _topol_struct &topol_struct,
                 double * box,
                 _coord_vel_force_struct &coord_vel_force_struct,
                 bool all) {

    int start;
    if (all) start = 1;
    else     start = topol_struct.Natm_biomol + 1;

    //Periodic boundary for crowders: this code performs similarly as iwrap in AMBER
    for (int k = start; k <= topol_struct.Natm; k++) {
        if (coord_vel_force_struct.coordx [k] < 0) {
            coord_vel_force_struct.coordx [k]     += box[0]; 
            coord_vel_force_struct.coord_oldx [k] += box[0];

        } else if (coord_vel_force_struct.coordx [k] > box[0]) {
            coord_vel_force_struct.coordx [k]     -= box[0];
            coord_vel_force_struct.coord_oldx [k] -= box[0];
        }

        if (coord_vel_force_struct.coordy [k] < 0) {
            coord_vel_force_struct.coordy [k]     += box[1]; 
            coord_vel_force_struct.coord_oldy [k] += box[1];

        } else if (coord_vel_force_struct.coordy [k] > box[1]) {
            coord_vel_force_struct.coordy [k]     -= box[1];
            coord_vel_force_struct.coord_oldy [k] -= box[1];
        }

        if (coord_vel_force_struct.coordz [k] < 0) {
            coord_vel_force_struct.coordz [k]     += box[2]; 
            coord_vel_force_struct.coord_oldz [k] += box[2];

        } else if (coord_vel_force_struct.coordz [k] > box[2]) {
            coord_vel_force_struct.coordz [k]     -= box[2];
            coord_vel_force_struct.coord_oldz [k] -= box[2];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void move_rigid_units (long int step,// FILE * f_error,
                       _simu_struct            &simu_struct,
                       const _topol_struct     &topol_struct,
                       _coord_vel_force_struct &coord_vel_force_struct) {

    for (int k = 1; k <= topol_struct.Natm; k++)
        leap_frog_move_atom (k, simu_struct.temp, simu_struct.dt, topol_struct, coord_vel_force_struct);

    if (step % simu_struct.step_adj_box == 0)
        adjust_box (topol_struct.bead_center, topol_struct.Natm, simu_struct.box, coord_vel_force_struct);

    wrap_coord (topol_struct, simu_struct.box, coord_vel_force_struct, 1);
}

///////////////////////////////////////////////////////////////////////////////////////////
void move_crowders (//long int step, FILE * f_error,
                    _simu_struct            &simu_struct,
                    const _topol_struct     &topol_struct,
                    _coord_vel_force_struct &coord_vel_force_struct) {

    for (int k = topol_struct.Natm_biomol + 1; k <= topol_struct.Natm; k++)
        leap_frog_move_atom (k, simu_struct.temp, simu_struct.dt, topol_struct, coord_vel_force_struct);

    wrap_coord (topol_struct, simu_struct.box, coord_vel_force_struct, 0);
}

#endif
