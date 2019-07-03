#ifndef EW_H
#define EW_H

void tabulate_CT (_EW_struct &EW_struct) {
    double C = 2.0 / sqrt (PI);
    for (int i = 0; i <= EW_struct.CTM; i++) {
        double j = 1e-3 * i;
        double term = erfc (j);
        EW_struct.CTE_array.push_back (term);
        EW_struct.CTF_array.push_back (term + C*j*exp (- j*j));
    }

/*    for (int i = 0; i <= 1571; i++) {
        double j = 1e-3 * i;
        EW_struct.sin_array.push_back (sin (j));
        EW_struct.cos_array.push_back (cos (j));
    }*/
}

///////////////////////////////////////////////////////////////
void prepare_fourier_space (const _simu_struct &simu_struct,
                            int Natm,
                            _EW_struct         &EW_struct) {
/*    int i, j, k, l;
    int S1, S2, S3;
    int km1, km2, km3, count;
    int np_busy, fOFFSET;
    int * index1 = NULL, * index2 = NULL, * index3 = NULL;
    double * kx_temp = NULL, * ky_temp = NULL, * kz_temp = NULL;
    double * const_temp = NULL;
    double p0, p1, p2, p3;
    double side_MIN, alpha;
    double kx, ky, kz, k2;*/
    //////////////////////////////////////////////////////////////////////////
/*    if (EW_struct.fLOAD) {
        for (int i = 1; i <= EW_struct.fLOAD; i++) {
            free (EW_struct.kosinus [i]);
            free (EW_struct.sinus   [i]);
        }
        free (EW_struct.kosinus); free (EW_struct.sinus);// free (A1); free (A2);
//        free (CONST_EWALD); free (fourier_x); free (fourier_y); free (fourier_z);
        kosinus = NULL; sinus = NULL;// A1 = NULL; A2 = NULL;
//        CONST_EWALD = NULL; fourier_x = NULL; fourier_y = NULL; fourier_z = NULL;
    }*/
    //////////////////////////////////////////////////////////////////////////
    double alpha = 1.0 / (EW_struct.sigma * EW_struct.sigma);
    //////////////////////////////////////////////////////////////////////////
    double tmp = EW_struct.kappa * sqrt (EW_struct.EXP_CUTOFF) / PI;
    int km1 = (int) (tmp * simu_struct.box[0]) + 1;
    int km2 = (int) (tmp * simu_struct.box[1]) + 1;
    int km3 = (int) (tmp * simu_struct.box[2]) + 1;
    //////////////////////////////////////////////////////////////////////////
    double p0 = 8.0 * PI * simu_struct.lB_T / (simu_struct.box[0] * simu_struct.box[1] * simu_struct.box[2]);
    double p1 = 2.0 * PI / simu_struct.box[0];
    double p2 = 2.0 * PI / simu_struct.box[1];
    double p3 = 2.0 * PI / simu_struct.box[2];
//    std::vector<double> const_temp, kx_temp, ky_temp, kz_temp;
    std::vector<int> index1, index2, index3;

/*    const_temp.push_back (0);
    kx_temp.push_back    (0);
    ky_temp.push_back    (0);
    kz_temp.push_back    (0);*/
    EW_struct.CONST_EWALD.push_back (0);
    EW_struct.fourier_x.push_back   (0);
    EW_struct.fourier_y.push_back   (0);
    EW_struct.fourier_z.push_back   (0);
    index1.push_back (0);
    index2.push_back (0);
    index3.push_back (0);

//    int count = 0;
    for (int i = -km1; i <= km1; i++)
        for (int j = -km2; j <= km2; j++)
            for (int k = -km3; k <= km3; k++) {
                if (i == 0 && j == 0 && k == 0) continue;
                int l = 0;
                for (l = 1; l <= EW_struct.N_vec; l++) {
                    int S1 = i + index1 [l];
                    int S2 = j + index2 [l];
                    int S3 = k + index3 [l];
                    if (S1 == 0 && S2 == 0 && S3 == 0) { l = 0; break; }
                }
                if (!l) continue;
                double kx = p1 * i;
                double ky = p2 * j;
                double kz = p3 * k;
                double k2 = kx*kx + ky*ky + kz*kz;
                double side_MIN = 0.25 * k2 * EW_struct.sigma * EW_struct.sigma;
                if (side_MIN < EW_struct.EXP_CUTOFF) {
//                    count++;
                    EW_struct.N_vec++;
//                    const_temp = (double *) realloc (const_temp, (count + 1) * sizeof (double));
//                    const_temp [count] = p0 * exp (-side_MIN) / k2;
//                    const_temp.push_back (p0 * exp (-side_MIN) / k2);
                    EW_struct.CONST_EWALD.push_back (p0 * exp (-side_MIN) / k2);

//                    index1 = (int *) realloc (index1, (count + 1) * sizeof (int));
//                    kx_temp = (double *) realloc (kx_temp, (count + 1) * sizeof (double));
//                    index1 [count] = i; kx_temp [count] = kx;
                    index1.push_back (i);
//                    kx_temp.push_back (kx);
                    EW_struct.fourier_x.push_back (kx);

//                    index2 = (int *) realloc (index2, (count + 1) * sizeof (int));
//                    ky_temp = (double *) realloc (ky_temp, (count + 1) * sizeof (double));
//                    index2 [count] = j; ky_temp [count] = ky;
                    index2.push_back (j);
//                    ky_temp.push_back (ky);
                    EW_struct.fourier_y.push_back (ky);

//                    index3 = (int *) realloc (index3, (count + 1) * sizeof (int));
//                    kz_temp = (double *) realloc (kz_temp, (count + 1) * sizeof (double));
//                    index3 [count] = k; kz_temp [count] = kz;
                    index3.push_back (k);
//                    kz_temp.push_back (kz);
                    EW_struct.fourier_z.push_back (kz);
                }
            }

//    free (index1); free (index2); free (index3);
    //////////////////////////////////////////////////////////////////////////
/*    EW_struct.N_vec = count;
    np_busy = N_vec % comm_size;
    if (my_rank < np_busy) {
        fLOAD = N_vec / comm_size + 1;
        fOFFSET = my_rank * fLOAD;
    } else {
        fLOAD = N_vec / comm_size;
        fOFFSET = my_rank * fLOAD + np_busy;
    }*/

//    kosinus = (double **) calloc (fLOAD + 1, sizeof (double *));
//    sinus   = (double **) calloc (fLOAD + 1, sizeof (double *));
    EW_struct.kosinus = (double **) calloc (EW_struct.N_vec + 1, sizeof (double *));
    EW_struct.sinus   = (double **) calloc (EW_struct.N_vec + 1, sizeof (double *));

    ///////// Need to loop through only charged beads? ////////
//    for (int i = 1; i <= fLOAD; i++) {
    printf ("Nvec  %d\n", EW_struct.N_vec);
    for (int i = 1; i <= EW_struct.N_vec; i++) {
        EW_struct.kosinus [i] = (double *) calloc (Natm + 1, sizeof (double));
        EW_struct.sinus   [i] = (double *) calloc (Natm + 1, sizeof (double));
    }

/*    A1 = (double *) calloc (fLOAD + 1, sizeof (double));
    A2 = (double *) calloc (fLOAD + 1, sizeof (double));
    CONST_EWALD = (double *) calloc (fLOAD + 1, sizeof (double));
    fourier_x   = (double *) calloc (fLOAD + 1, sizeof (double));
    fourier_y   = (double *) calloc (fLOAD + 1, sizeof (double));
    fourier_z   = (double *) calloc (fLOAD + 1, sizeof (double));*/

    EW_struct.A1.resize (EW_struct.N_vec + 1, 0);
    EW_struct.A2.resize (EW_struct.N_vec + 1, 0);

    //////////////////////////////////////////////////////////////////////////
/*    for (int i = 1; i <= fLOAD; i++) {
        CONST_EWALD [i] = (const_temp + fOFFSET) [i];
        fourier_x   [i] = (kx_temp    + fOFFSET) [i];
        fourier_y   [i] = (ky_temp    + fOFFSET) [i];
        fourier_z   [i] = (kz_temp    + fOFFSET) [i];
    }*/
//    free (kx_temp); free (ky_temp); free (kz_temp);
    //////////////////////////////////////////////////////////////////////////
    double total = 0;
    for (int i = 1; i <= EW_struct.N_vec; i++)
//        total += EW_struct.const_temp [i];
        total += EW_struct.CONST_EWALD [i];
    EW_struct.CONST_ADD    =  0.5*total - simu_struct.lB_T * sqrt (alpha / PI);
    EW_struct.CONST_REMOVE = -0.5*total - simu_struct.lB_T * sqrt (alpha / PI);
//    free (const_temp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepare_kosinus_sinus (const _topol_struct           &topol_struct,
                            _EW_struct                    &EW_struct,
                            const _coord_vel_force_struct &coord_vel_force_struct) {
//    for (int i = 1; i <= fLOAD; i++) {
    for (int i = 1; i <= EW_struct.N_vec; i++) {
        EW_struct.A1 [i] = 0;    EW_struct.A2 [i] = 0;
        for (int j = 1; j <= topol_struct.Natm; j++) {
            int J1 = topol_struct.maxi_key [j];
            if (topol_struct.QC [J1]) {
//                double tmp = kx * x[j] + ky * y[j] + kz * z[j];
                double tmp = EW_struct.fourier_x [i] * coord_vel_force_struct.coordx [j] + EW_struct.fourier_y [i] * coord_vel_force_struct.coordy [j] + EW_struct.fourier_z [i] * coord_vel_force_struct.coordz [j];
                ///////////////////////////////////////////////////////
                EW_struct.kosinus [i][j] = cos (tmp);
                EW_struct.sinus   [i][j] = sin (tmp);
                ///////////////////////////////////////////////////////
                EW_struct.A1 [i] += topol_struct.QC [J1] * EW_struct.kosinus [i][j];
                EW_struct.A2 [i] += topol_struct.QC [J1] * EW_struct.sinus   [i][j];
            }
        }
    }
}

/////////////////////////////////////////////////
void Coulomb_interactions_fourier (const _topol_struct     &topol_struct,
                                   _coord_vel_force_struct &coord_vel_force_struct,
                                   const _EW_struct        &EW_struct) {
    for (int i = 1; i <= EW_struct.N_vec; i++)
        for (int j = 1; j <= topol_struct.Natm; j++) {
            int J1 = topol_struct.maxi_key [j];
            if (topol_struct.QC [J1]) {
                double tmp = topol_struct.QC [J1] * (EW_struct.A1 [i] * EW_struct.sinus [i][j] - EW_struct.A2 [i] * EW_struct.kosinus [i][j]);
                coord_vel_force_struct.flx [j] += EW_struct.CONST_EWALD [i] * EW_struct.fourier_x [i] * tmp;
                coord_vel_force_struct.fly [j] += EW_struct.CONST_EWALD [i] * EW_struct.fourier_y [i] * tmp;
                coord_vel_force_struct.flz [j] += EW_struct.CONST_EWALD [i] * EW_struct.fourier_z [i] * tmp;
            }
        }
}

//////////////////////////////////////////////////////////////
/*double COULOMB_TRUNCATED_ENERGY_EWALD (double r2, int Jn12) {
    if (r2 > D2_CT_CUTOFF) return (0);
    else {
        double r1, A;
        r1 = sqrt (r2);
        A = Q_C [Jn12] * lB / r1;
        r1 = kappa * r1;
        return (A * CTE (r1)); 
    }
}

///////////////////////////////////////////////////////////////////////
double COULOMB_TRUNCATED_FORCE_EWALD (double &r2, double Q_C,
                                      const _simu_struct &simu_struct,
                                      const _EW_struct   &EW_struct) {
    if (r2 > simu_struct.D2_CT_CUTOFF) return (0);
    else {
        double r1 = sqrt (r2); 
        double A = Q_C * simu_struct.lB_T;// / r1;
        r1 *= EW_struct.kappa;
        return A * EW_struct.CTF_array [int (1000*r1)] / r2;
    }
}*/

#endif
