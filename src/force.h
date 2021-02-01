#ifndef FORCE_H
#define FORCE_H

/////////////////////////////////////////////////////////////////////////////////
/*void forces_PRMD (_simu_struct            &simu_struct,
                  _topol_struct           &topol_struct,
                  _coord_vel_force_struct &coord_vel_force_struct,
                  const _EW_struct        &EW_struct) {

    if (simu_struct.relist_step % simu_struct.LR_skip == 0) {
        for (int i = 1; i <= topol_struct.Natm; i++) {
            coord_vel_force_struct.flx [i] = 0;
            coord_vel_force_struct.fly [i] = 0;
            coord_vel_force_struct.flz [i] = 0;
        }
        LR_interactions_lists (simu_struct, topol_struct, coord_vel_force_struct, EW_struct);
//        Coulomb_interactions_fourier (topol_struct, coord_vel_force_struct, EW_struct);
    }

    for (int i = 1; i <= topol_struct.Natm; i++) {
        coord_vel_force_struct.forcex [i] = coord_vel_force_struct.flx [i];
        coord_vel_force_struct.forcey [i] = coord_vel_force_struct.fly [i];
        coord_vel_force_struct.forcez [i] = coord_vel_force_struct.flz [i];
    }

    SR_interactions_lists (simu_struct, topol_struct, coord_vel_force_struct, EW_struct);
}*/

///////////////////////////////////////////////////
//// Compute force from PMF for explicit ions
///////////////////////////////////////////////
/*double compute_force_qq (double **pair_potential,
                         double r, int i1, int i2) {
    int index = 0;
    if      (i1 + i2 == 17)            index = 2;
    else if (i1 + i2 == 10)            index = 4;
    else if (i1 + i2 == 15)            index = 6;
    else if (i1 + i2 == 8)             index = 7;
    else if (i1 + i2 == 18)            index = 1;
    else if (i1 == 8 && i2 == 8)       index = 8;
    else if (i1 + i2 == 9)             index = 9;
    else if (i1 + i2 == 16 && i1 != 8) index = 3;
    else if (i1 == 7 && i2 == 7)       index = 5;

    double d = 40 * (r - pair_potential [index][0]); // 40 = 1/0.025
    int i = std::floor (d) + 1;
    if (i > 1)
        return 40 * (pair_potential [index][i] - pair_potential [index][i+1]); // 40 = 1/0.025
    else {
        std::cout << "There are two particles approaching too close !!!  " << i1 << " " << i2 << std::endl;
        std::cout << "Distance " << r << std::endl;
        exit (0);
    }
}*/

///////////////////////////////////////////////////////////////////////////
//  Update new Ze to take into account explicit Mg
//  V assumes infinitely long charge, from Wilson et al BJ, 1980, 30, 317
//////////////////////////////////////////////////////////////////////////
void update_Ze (_simu_struct &simu_struct,
                double       *ion_conc) {

    double b3 = simu_struct.b * simu_struct.b * simu_struct.b;

    /// V assumes infinitely long charge, taken from Wilson et al BJ, 1980, 30, 317
    double V2 = 4*PI*E * b3 * (1/simu_struct.Ze - 1./simu_struct.ion_q) * (simu_struct.ion_q + 1);
    double V1 = 4*PI*E * b3 * (1/simu_struct.Ze - 1) * 2;
    double C1V1 = ion_conc [2] * V1 * 6.022e-04;
    double C2V2 = ion_conc [0] * V2 * 6.022e-04;

    // Assuming total condensation contribution from ions equal to (1 - Ze)
    double theta1 = 0;
    if (simu_struct.ion_q == 2)
        theta1 = (sqrt( pow(C1V1, 4) + 8*E * C1V1*C1V1 * C2V2 * (1. - simu_struct.Ze)) - C1V1*C1V1) / (4*E * C2V2);
    else if (simu_struct.ion_q == 3) {
        double C1V1_3 = C1V1*C1V1*C1V1;
        double C1V1_6 = C1V1_3 * C1V1_3;
        double C1V1_9 = C1V1_3 * C1V1_6;

        double C2V2_2 = C2V2*C2V2;
        double C2V2_3 = C2V2_2 * C2V2;

        double theta_sqr = (1. - simu_struct.Ze) * (1. - simu_struct.Ze);

        double tmp = sqrt (4* C1V1_9 * C2V2_3 + 81*E*E* C1V1_6 * C2V2_2*C2V2_2 * theta_sqr);
        tmp += 9*E* C1V1_3 * C2V2_2 * (1. - simu_struct.Ze);
        tmp = cbrt(tmp);
        theta1 = (tmp / (1.26*E*C2V2) - 1.26* C1V1_3 / (E*tmp)) / 3; /// 1.26 = cbrt(2)
    } else {
        printf ("Ion charge should be either +2 or +3!!!!!\n");
        exit (0);
    }

    simu_struct.Ze = 1. - theta1;
}

/////////////////////////////////////////////////////////////////////////////////
///// Merge the long-range part of the pair potential into the Debye-Huckel curve
///////////////////////////////////////////////////////////////////////////////
void update_pmf_ion_P_implicit (std::vector<double> &pair_potential,
                                const _simu_struct &simu_struct) {
    double r = pair_potential [0];
    int i = 1;
    while (true) {
        double DH = -simu_struct.ion_q * simu_struct.lB_T * simu_struct.Ze * exp (- r*simu_struct.kappaD) / r;
        double error = 0;
        if      (simu_struct.ion_q == 2) error = (DH - pair_potential [i]) * exp (-25. / (r*r));
        else if (simu_struct.ion_q == 3) error = (DH - pair_potential [i]) * exp (-150. / (r*r*r));
        else {
            printf ("Ion charge should be either +2 or +3 !!!\n");
            exit (0);
        }
        //if (r < 6.) printf (" r = %8.3f    old-pot = %12.3e", r, pair_potential [i]);
        pair_potential [i] += error;
        //if (r < 6.) printf ("  new-pot = %12.3e\n", pair_potential [i]);
        r += 0.025;
        if (r > 40.) break;
        i++;
    }
}

///////////////////////////////////////////////////////////////////////////
double compute_force_ion_P (const std::vector<double> &pair_potential,
                            double r, double &E_Q) {
    double d = 40 * (r - pair_potential [0]); // 40 = 1/0.025
    int i = std::ceil (d);
    if (i > 1) {
        double r1 = pair_potential [0] + (i - 1)*0.025;
        double slope = 40 * (pair_potential [i+1] - pair_potential [i]);  // 40 = 1/0.025
        E_Q += pair_potential [i] + slope * (r - r1);
        return -slope;
    } else {
        std::cout << "There are two particles approaching too close !!! Distance " << r << std::endl;
        exit (0);
    }
}

/////////////////////////////////////////////////////////////////////////
/*void LR_interactions_lists (_simu_struct            &simu_struct,
                            const _topol_struct     &topol_struct,
                            _coord_vel_force_struct &coord_vel_force_struct,
                            const _EW_struct        &EW_struct) {

    for (int k = 1; k <= simu_struct.list_mass_LR; k += 2) {
        int i = simu_struct.list_content_LR [k];
        int j = simu_struct.list_content_LR [k + 1];

        int i1 = topol_struct.maxi_key [i];
        int i2 = topol_struct.maxi_key [j];

        int j1;
        if (i1 < i2)
             j1 = (i1 - 1) * NA - (i1 - 1) * i1 / 2 + i2;
        else j1 = (i2 - 1) * NA - (i2 - 1) * i2 / 2 + i1;

        double vec[3];
        vec[0] = coord_vel_force_struct.coordx [i] - coord_vel_force_struct.coordx [j]; 
        vec[1] = coord_vel_force_struct.coordy [i] - coord_vel_force_struct.coordy [j]; 
        vec[2] = coord_vel_force_struct.coordz [i] - coord_vel_force_struct.coordz [j]; 
        half_shift (vec, simu_struct.box);

        double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];

//        double tmp = COULOMB_TRUNCATED_FORCE (r2, topol_struct.Q_C [j1], simu_struct.lB_T, simu_struct.D2_CT_CUTOFF);
//        double tmp = COULOMB_TRUNCATED_FORCE_EWALD (r2, topol_struct.Q_C [j1], simu_struct, EW_struct);
        double tmp = COULOMB_TRUNCATED_FORCE_WF (r2, topol_struct.Q_C [j1], simu_struct.lB_T, simu_struct.D2_CT_CUTOFF);

        vec[0] = vec[0] * tmp;
        vec[1] = vec[1] * tmp;
        vec[2] = vec[2] * tmp;

        coord_vel_force_struct.flx [i] += vec[0];
        coord_vel_force_struct.fly [i] += vec[1];
        coord_vel_force_struct.flz [i] += vec[2];

        coord_vel_force_struct.flx [j] -= vec[0];
        coord_vel_force_struct.fly [j] -= vec[1];
        coord_vel_force_struct.flz [j] -= vec[2];
    }
}

////////////////////////////////////////////////////////////////////////
void SR_interactions_lists (_simu_struct            &simu_struct,
                            _topol_struct           &topol_struct,
                            _coord_vel_force_struct &coord_vel_force_struct,
                            const _EW_struct        &EW_struct) {

    for (int k = 1; k <= simu_struct.list_mass_SR; k += 2) {
        int i = simu_struct.list_content_SR [k];
        int j = simu_struct.list_content_SR [k + 1];

        int i1 = topol_struct.maxi_key [i];
        int i2 = topol_struct.maxi_key [j];

        int j1;
        if (i1 < i2)
             j1 = (i1 - 1) * NA - (i1 - 1) * i1 / 2 + i2;
        else j1 = (i2 - 1) * NA - (i2 - 1) * i2 / 2 + i1;

        double vec [3];
        vec [0] = coord_vel_force_struct.coordx [i] - coord_vel_force_struct.coordx [j];
        vec [1] = coord_vel_force_struct.coordy [i] - coord_vel_force_struct.coordy [j];
        vec [2] = coord_vel_force_struct.coordz [i] - coord_vel_force_struct.coordz [j];
        half_shift (vec, simu_struct.box);

        double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
        double tmp = LJ_FORCE (r2, topol_struct.D_LJ [j1], topol_struct.F_LJ [j1], topol_struct.D2_LJ_CUTOFF [j1]);

        i1 = topol_struct.part_key [i];
        i2 = topol_struct.part_key [j];

        if (i1 != 1 && i1 != 2 && i2 != 1 && i2 != 2)
//            tmp += COULOMB_TRUNCATED_FORCE (r2, topol_struct.Q_C [j1], simu_struct.lB_T, simu_struct.D2_CT_CUTOFF);
//            tmp += COULOMB_TRUNCATED_FORCE_EWALD (r2, topol_struct.Q_C [j1], simu_struct, EW_struct);
            tmp += COULOMB_TRUNCATED_FORCE_WF (r2, topol_struct.Q_C [j1], simu_struct.lB_T, simu_struct.D2_CT_CUTOFF);

        vec[0] = vec[0] * tmp;
        vec[1] = vec[1] * tmp;
        vec[2] = vec[2] * tmp;

        coord_vel_force_struct.forcex [i] += vec[0];
        coord_vel_force_struct.forcey [i] += vec[1];
        coord_vel_force_struct.forcez [i] += vec[2];

        coord_vel_force_struct.forcex [j] -= vec[0];
        coord_vel_force_struct.forcey [j] -= vec[1];
        coord_vel_force_struct.forcez [j] -= vec[2];
    }
}*/

/////////////////////////////////////////////////////////////////////////////
double COULOMB_TRUNCATED_FORCE (double &r2, double Q_C,
                                double &lB, double &D2_CUTOFF) {
    if (r2 > D2_CUTOFF) return (0);
    else return Q_C * lB / r2;
/*    else {
        double r3 = r2 * sqrt(r2);
        return (Q_C * lB / r3);
    }*/
}

//////////////////////////////////////////////////////////////////////
// Use truncation-based to deal with long-ranged electrostatics
// This approach is the so-called shift-potential or shift-force
// (depending on whether the potential or force at the cutoff distance is set to 0).
// Here we use shift force, which was shown to reproduce EW force very well
// See more at Wolf et al, J Chem Phys 1999, 110, 8254
//             Fennell and Gezelter, J Chem Phys 2006, 124, 234104
//             Hansen et al, J Phys Chem B 2012, 116, 5738
// Note that CHARMM also has similar feature, but the functional form is different
// and thus the force is not accurate as this approach here
///////////////////////////////////////////////////////////////////
double COULOMB_TRUNCATED_FORCE_WF (double &r2, double Q_C,
                                   double &lB, double &D2_CUTOFF) {
    if (r2 >= D2_CUTOFF) return (0);
    else {
        double invR = 1./r2 - 1./D2_CUTOFF;
//        double r1 = sqrt (r2);
//        double R1 = sqrt (D2_CT_CUTOFF);
//        double dr = R1 - r1;
//        return (Q_C * lB / r1 * (invR - 2./(r2*r1)*dr));
        return (Q_C * lB * invR);
    }
}

///////////////////////////////////////////////////////////////////////
double DEBYE_HUCKEL_FORCE (double &r, double &kappaD,
                           double E_DH, double &D_CUTOFF,
                           double &E_Q) {
    if (r >= D_CUTOFF) return (0);
    else {
    //    r3 = r2 * r1;
        double x = kappaD * r;
        double E0 = E_DH * exp (- x);
        E_Q += E0 / r;
        return E0 * (1 + x) / (r*r);
    }
}

///////////////////////////////////////////////////////////////////////
double LJ_FORCE (double &x, double D_LJ, double F_LJ, double &E_EXCLUDED) {
    if (x > D_LJ) return (0);
    else {
//        double r1 = sqrt (x2);
        double dr = x + 1.5852 - D_LJ;  // Check HARD-CODE
/*        double r4 = 1.5852 / dr;
        double r2 = r4 * r4;
        r4 = r2 * r2;
        double r8 = r4 * r4;
        double r14 = r8 * r4 * r2;

        return F_LJ * (r14 - r8) * dr / r1;*/

        double r2 = 1.5852 / dr;
        r2 = r2 * r2;
        double r6 = r2 * r2 * r2;
        double r12 = r6 * r6;

        E_EXCLUDED += 0.0833333333 * F_LJ * (r12 - 2*r6 + 1.);  /// 1/12 * F_LJ = E_LJ
        return F_LJ * (r12 - r6) / dr;
    }
}

/////////////////////////////////////////////////////////////////////////////
void non_bonded_interaction (_simu_struct            &simu_struct,
                             const _topol_struct     &topol_struct,
                             _coord_vel_force_struct &coord_vel_force_struct,
                             _energy_struct          &energy_struct) {
    /// Electrostatic
    for (int k = 1; k <= simu_struct.DHneighborlist_mass; k += 2) {
        int i = simu_struct.DHneighborlist [k];
        int j = simu_struct.DHneighborlist [k + 1];

        int i1 = topol_struct.maxi_key [i];
        int i2 = topol_struct.maxi_key [j];

        double vec [3];
        vec[0] = coord_vel_force_struct.coordx [i] - coord_vel_force_struct.coordx [j];
        vec[1] = coord_vel_force_struct.coordy [i] - coord_vel_force_struct.coordy [j];
        vec[2] = coord_vel_force_struct.coordz [i] - coord_vel_force_struct.coordz [j];
        half_shift (vec, simu_struct.box);

        double r = sqrt (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
        
        if (r > simu_struct.D_CUTOFF) continue;

        double tmp = 0;

        if (((i1 == 1 && i2 >= 7) || (i1 >= 7 && i2 == 1)) && (r < 40.))
            tmp = compute_force_ion_P (simu_struct.pair_potential, r, energy_struct.E_Q);
        else {
            //if (i1 == 1 && i2 == 1) printf ("P-P force\n");
            int j1;
            if (i1 < i2)
                 j1 = (i1 - 1) * NA - (i1 - 1) * i1 / 2 + i2;
            else j1 = (i2 - 1) * NA - (i2 - 1) * i2 / 2 + i1;

            tmp = DEBYE_HUCKEL_FORCE (r, simu_struct.kappaD, topol_struct.E_DH [j1], simu_struct.D_CUTOFF, energy_struct.E_Q);
        }

        tmp /= r;

        vec[0] = vec[0] * tmp;
        vec[1] = vec[1] * tmp;
        vec[2] = vec[2] * tmp;

        coord_vel_force_struct.forcex [i] += vec[0];
        coord_vel_force_struct.forcey [i] += vec[1];
        coord_vel_force_struct.forcez [i] += vec[2];

        coord_vel_force_struct.forcex [j] -= vec[0];
        coord_vel_force_struct.forcey [j] -= vec[1];
        coord_vel_force_struct.forcez [j] -= vec[2];
    }

    // LJ
    for (int k = 1; k <= simu_struct.HBneighborlist_mass; k += 2) {
        int i = simu_struct.HBneighborlist [k];
        int j = simu_struct.HBneighborlist [k + 1];

        int i1 = topol_struct.maxi_key [i];
        int i2 = topol_struct.maxi_key [j];

        int j1;
        if (i1 < i2)
             j1 = (i1 - 1) * NA - (i1 - 1) * i1 / 2 + i2;
        else j1 = (i2 - 1) * NA - (i2 - 1) * i2 / 2 + i1;

        double vec [3];
        vec[0] = coord_vel_force_struct.coordx [i] - coord_vel_force_struct.coordx [j];
        vec[1] = coord_vel_force_struct.coordy [i] - coord_vel_force_struct.coordy [j];
        vec[2] = coord_vel_force_struct.coordz [i] - coord_vel_force_struct.coordz [j];
        half_shift (vec, simu_struct.box);

        double r = sqrt (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
        double tmp = LJ_FORCE (r, topol_struct.D_LJ [j1], topol_struct.F_LJ [j1], energy_struct.E_EXCLUDED);
        tmp /= r;

        vec[0] = vec[0] * tmp;
        vec[1] = vec[1] * tmp;
        vec[2] = vec[2] * tmp;

        coord_vel_force_struct.forcex [i] += vec[0];
        coord_vel_force_struct.forcey [i] += vec[1];
        coord_vel_force_struct.forcez [i] += vec[2];

        coord_vel_force_struct.forcex [j] -= vec[0];
        coord_vel_force_struct.forcey [j] -= vec[1];
        coord_vel_force_struct.forcez [j] -= vec[2];
    }

  //  for (long int k = 1; k <= simu_struct.neighborlist_mass; k += 2) {
  //      int i = simu_struct.neighborlist [k];
  //      int j = simu_struct.neighborlist [k + 1];

  //      int i1 = topol_struct.maxi_key [i];
  //      int i2 = topol_struct.maxi_key [j];

  //      int j1;
  //      if (i1 < i2)
  //           j1 = (i1 - 1) * NA - (i1 - 1) * i1 / 2 + i2;
  //      else j1 = (i2 - 1) * NA - (i2 - 1) * i2 / 2 + i1;

  //      double vec[3];
  //      vec[0] = coord_vel_force_struct.coordx [i] - coord_vel_force_struct.coordx [j];
  //      vec[1] = coord_vel_force_struct.coordy [i] - coord_vel_force_struct.coordy [j];
  //      vec[2] = coord_vel_force_struct.coordz [i] - coord_vel_force_struct.coordz [j];
  //      half_shift (vec, simu_struct.box);

  //      double r = sqrt (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  //      double tmp = 0;

  //      //// Ion - P case
  //      if ((i1 == 1 && i2 >= 7) || (i1 >= 7 && i2 == 1)) {
  //          if (not std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), i) && \
  //              not std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), j)) {
  //              if (r < 40.)
  //                  tmp += compute_force_ion_P (simu_struct.pair_potential, r, energy_struct.E_Q);
  //              else
  //                  tmp += DEBYE_HUCKEL_FORCE (r, simu_struct.kappaD, topol_struct.E_DH [j1], simu_struct.D_CUTOFF, energy_struct.E_Q);
  //          } else tmp += LJ_FORCE (r, topol_struct.D_LJ [j1], topol_struct.F_LJ [j1], energy_struct.E_EXCLUDED);

  //      /// Other cases
  //      } else {
  //          tmp = LJ_FORCE (r, topol_struct.D_LJ [j1], topol_struct.F_LJ [j1], energy_struct.E_EXCLUDED);
  //          if ((topol_struct.Q_C [j1] > 1e-3 || topol_struct.Q_C [j1] < -1e-3) && \
  //               not std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), i) && \
  //               not std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), j))
  //             /* if (i1 + i2 == 8 && r < 40.)
  //                  tmp += compute_force_Mg_P (simu_struct.pair_potential [7], r, energy_struct.E_Q);
  //              else*/
  //              tmp += DEBYE_HUCKEL_FORCE (r, simu_struct.kappaD, topol_struct.E_DH [j1], simu_struct.D_CUTOFF, energy_struct.E_Q);
  //      }
/*//        if (simu_struct.Coulomb) {
  //          if (r < topol_struct.D_LJ [j1])
  //              tmp += LJ_FORCE (r, topol_struct.D_LJ [j1], topol_struct.F_LJ [j1]);
  //          if (topol_struct.Q_C [j1] != 0)
  //              tmp += COULOMB_TRUNCATED_FORCE (r2, topol_struct.Q_C [j1], simu_struct.lB_T, simu_struct.D2_CUTOFF);
  //      } else {
  //          // Use numerical forces
  //          if (topol_struct.Q_C [j1] != 0)
  //              if ((i1 != 1 || i2 != 1) && r < 50) {
  //                  tmp = compute_force_qq (simu_struct.pair_potential, r, i1, i2);
  //                  tmp -= topol_struct.Q_C [j1] * simu_struct.lB_T_D2;
  //                  // Normal Coulomb shifted_force for P-P
  //              } else tmp = COULOMB_TRUNCATED_FORCE_WF (r2, topol_struct.Q_C [j1], simu_struct.lB_T, simu_struct.D2_CUTOFF);

  //          // LJ for P-P and non-charged pairs
  //          if (topol_struct.Q_C [j1] == 0 || (i1 == 1 && i2 == 1))
  //              tmp += LJ_FORCE (r, topol_struct.D_LJ [j1], topol_struct.F_LJ [j1]);
  //      }*/

  //      tmp /= r;

  //      vec[0] = vec[0] * tmp;
  //      vec[1] = vec[1] * tmp;
  //      vec[2] = vec[2] * tmp;

  //      coord_vel_force_struct.forcex [i] += vec[0];
  //      coord_vel_force_struct.forcey [i] += vec[1];
  //      coord_vel_force_struct.forcez [i] += vec[2];

  //      coord_vel_force_struct.forcex [j] -= vec[0];
  //      coord_vel_force_struct.forcey [j] -= vec[1];
  //      coord_vel_force_struct.forcez [j] -= vec[2];
  //  }
}

/////////////////////////////////////////////////////////////////////
void gen_bond_force (int i1, int i2, double k0, double r0,
                     _coord_vel_force_struct &coord_vel_force_struct,
                     double    &E_BOND) {
    double v [3];
    v [0] = coord_vel_force_struct.coordx [i2] - coord_vel_force_struct.coordx [i1];
    v [1] = coord_vel_force_struct.coordy [i2] - coord_vel_force_struct.coordy [i1];
    v [2] = coord_vel_force_struct.coordz [i2] - coord_vel_force_struct.coordz [i1];

    double r = sqrt (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    double dr = r0 - r;
//    double tmp = k0 * dr / r;
//    E_BOND += tmp * dr;
    E_BOND += k0 * dr * dr;
    double tmp = 2.*k0 * dr / r;
//    tmp *= 2.;

    r = tmp * v [0];
    coord_vel_force_struct.forcex [i1] -= r;
    coord_vel_force_struct.forcex [i2] += r;

    r = tmp * v [1];
    coord_vel_force_struct.forcey [i1] -= r;
    coord_vel_force_struct.forcey [i2] += r;

    r = tmp * v [2];
    coord_vel_force_struct.forcez [i1] -= r;
    coord_vel_force_struct.forcez [i2] += r;
}

////////////////////////////////////////////////////////////////////////
void gen_valence_force (int i1, int i2, int i3, double k0, double theta0,
                        _coord_vel_force_struct &coord_vel_force_struct,
                        double  &E_ANGLE) {
    double v1 [3], v2 [3];
    int I [4];

    I [1] = i1; I [2] = i2; I [3] = i3;

    v1 [0] = coord_vel_force_struct.coordx [i1] - coord_vel_force_struct.coordx [i2];
    v1 [1] = coord_vel_force_struct.coordy [i1] - coord_vel_force_struct.coordy [i2];
    v1 [2] = coord_vel_force_struct.coordz [i1] - coord_vel_force_struct.coordz [i2];

    v2 [0] = coord_vel_force_struct.coordx [i3] - coord_vel_force_struct.coordx [i2];
    v2 [1] = coord_vel_force_struct.coordy [i3] - coord_vel_force_struct.coordy [i2];
    v2 [2] = coord_vel_force_struct.coordz [i3] - coord_vel_force_struct.coordz [i2];

    double v1_norm = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    v1_norm = 1. / v1_norm;

    double v2_norm = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
    v2_norm = 1. / v2_norm;

    double v12_norm = sqrt (v1_norm * v2_norm);

    double kosinus = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) * v12_norm;
    double sinus, theta;
    if (kosinus < -0.99995) {
        sinus = 0;
        theta = PI;
    } else {
        sinus = sqrt (1. - kosinus * kosinus);
        theta = acos (kosinus);
    }

    v1_norm = kosinus * v1_norm;
    v2_norm = kosinus * v2_norm;

    if (sinus) {
        double dtheta = theta0 - theta;
//        double tmp = k0 * dtheta / sinus;
        sinus = - 2.*k0 * dtheta / sinus;
//        sinus = -2.*tmp;
//        E_ANGLE += tmp * dtheta;
        E_ANGLE += k0 * dtheta * dtheta;

        for (int k = 0; k <= 2; k++)
            for (int l = 1; l <= 3; l++) {
                double tmp1 = Kdelta (l, 1);
                double tmp2 = Kdelta (l, 2);
                double tmp3 = Kdelta (l, 3);

                double C1 = tmp1 - tmp2;
                double C2 = tmp3 - tmp2;

                tmp1 = C1 * v2 [k] + C2 * v1 [k];
                tmp2 = C1 * v1 [k];
                tmp3 = C2 * v2 [k];

                C1 = tmp1 * v12_norm - tmp2 * v1_norm - tmp3 * v2_norm;

                tmp1 = sinus * C1;
                if (k == 0) coord_vel_force_struct.forcex [I [l]] += tmp1;
                if (k == 1) coord_vel_force_struct.forcey [I [l]] += tmp1;
                if (k == 2) coord_vel_force_struct.forcez [I [l]] += tmp1;
            }
    }
}

////////////////////////////////////////////////////////////////
void reduce_forces (int i1, int i2, double scale_factor,
                    double * box,
                    const _topol_struct     &topol_struct,
                    _coord_vel_force_struct &coord_vel_force_struct) {
    int Jn1 = topol_struct.maxi_key [i1];
    int Jn2 = topol_struct.maxi_key [i2];

    int Jn12;
    if (Jn1 < Jn2)
         Jn12 = (Jn1 - 1) * NA - (Jn1 - 1) * Jn1 / 2 + Jn2;
    else Jn12 = (Jn2 - 1) * NA - (Jn2 - 1) * Jn2 / 2 + Jn1;

    // CC_force_single_pair ( i1, i2, Jn12, scale_factor, LJ_FORCE );
    ///// Paste codes from CC_force_single_pair
    double vec [3];
    vec [0] = coord_vel_force_struct.coordx [i1] - coord_vel_force_struct.coordx [i2];
    vec [1] = coord_vel_force_struct.coordy [i1] - coord_vel_force_struct.coordy [i2];
    vec [2] = coord_vel_force_struct.coordz [i1] - coord_vel_force_struct.coordz [i2];

    half_shift (vec, box);

    double r = sqrt (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    double x = 0; // trash variable
    double tmp = scale_factor / r * LJ_FORCE (r, topol_struct.D_LJ [Jn12], topol_struct.F_LJ [Jn12], x); // Jn12 );

    vec[0] = vec[0] * tmp;
    vec[1] = vec[1] * tmp;
    vec[2] = vec[2] * tmp;

    coord_vel_force_struct.forcex [i1] += vec[0];
    coord_vel_force_struct.forcey [i1] += vec[1];
    coord_vel_force_struct.forcez [i1] += vec[2];

    coord_vel_force_struct.forcex [i2] -= vec[0];
    coord_vel_force_struct.forcey [i2] -= vec[1];
    coord_vel_force_struct.forcez [i2] -= vec[2];
}

//////////////////////////////////////////////////////////////////
/*void bead_bond_force (int i1, int i2, double k0, double r0,
                      double * box,
                      const _topol_struct     &topol_struct,
                      _coord_vel_force_struct &coord_vel_force_struct,
                      double &E_BOND) {
    gen_bond_force (i1, i2, k0, r0, coord_vel_force_struct, E_BOND);
    //reduce_forces (i1, i2, -1.0, box, topol_struct, coord_vel_force_struct);
}

///////////////////////////////////////////////////////////////////////////////
void bead_valence_force (int i1, int i2, int i3, double k0, double theta0,
                         double * box,
                         const _topol_struct     &topol_struct,
                         _coord_vel_force_struct &coord_vel_force_struct,
                         double &E_ANGLE) {
    gen_valence_force (i1, i2, i3, k0, theta0, coord_vel_force_struct, E_ANGLE);
    //reduce_forces (i1, i3, -1.0, box, topol_struct, coord_vel_force_struct);
}*/

////////////////////////////////////////////////////////////////////
void stacking_dihedral_force (int i1, int i2, int i3, int i4, double k0, double psi0,
                              _coord_vel_force_struct &coord_vel_force_struct) {

    int I [5];
    I [1] = i1;    I [2] = i2;       I [3] = i3;       I [4] = i4;

    double v1 [3], v2 [3], v3 [3];
    v1 [0] = coord_vel_force_struct.coordx [i2] - coord_vel_force_struct.coordx [i1];
    v1 [1] = coord_vel_force_struct.coordy [i2] - coord_vel_force_struct.coordy [i1];
    v1 [2] = coord_vel_force_struct.coordz [i2] - coord_vel_force_struct.coordz [i1];

    v2 [0] = coord_vel_force_struct.coordx [i3] - coord_vel_force_struct.coordx [i2];
    v2 [1] = coord_vel_force_struct.coordy [i3] - coord_vel_force_struct.coordy [i2];
    v2 [2] = coord_vel_force_struct.coordz [i3] - coord_vel_force_struct.coordz [i2];

    v3 [0] = coord_vel_force_struct.coordx [i4] - coord_vel_force_struct.coordx [i3];
    v3 [1] = coord_vel_force_struct.coordy [i4] - coord_vel_force_struct.coordy [i3];
    v3 [2] = coord_vel_force_struct.coordz [i4] - coord_vel_force_struct.coordz [i3];

    double m [3], n [3];
    cross_product2 (v1, v2, m);
    cross_product2 (v2, v3, n);

    double m_norm = m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
    m_norm = 1.0 / m_norm;

    double n_norm = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
    n_norm = 1.0 / n_norm;

    double mn_norm = sqrt (m_norm * n_norm);

    double tmp1    = sqrt (v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    double kosinus = (m[0]*n[0] + m[1]*n[1] + m[2]*n[2]) * mn_norm;
    double sinus   = (v1[0]*n[0] + v1[1]*n[1] + v1[2]*n[2]) * tmp1 * mn_norm;
    double psi     = atan2 (sinus, kosinus);

    m_norm = kosinus * m_norm;
    n_norm = kosinus * n_norm;

    sinus = - 2.0 * k0 * (psi0 - psi) / sinus;

    for (int k = 0; k <= 2; k++) {
        double A11 = 0, A22 = 0, A33 = 0, A12 = 0, A23 = 0, A13 = 0;

        for (int l = 0; l <= 2; l++) {
            double tmp1 = 1 - Kdelta (l, k);
            A11 += tmp1 * v1[l] * v1[l];
            A22 += tmp1 * v2[l] * v2[l];
            A33 += tmp1 * v3[l] * v3[l];
            A12 += tmp1 * v1[l] * v2[l];
            A23 += tmp1 * v2[l] * v3[l];
            A13 += tmp1 * v1[l] * v3[l];
        }

        for (int l = 1; l <= 4; l++) {
            double tmp1 = Kdelta (l, 1);
            double tmp2 = Kdelta (l, 2);
            double tmp3 = Kdelta (l, 3);
            double tmp4 = Kdelta (l, 4);

            double C1 = A22 * (tmp3 - tmp4) + A23 * (tmp3 - tmp2);
            double C2 = A12 * (tmp3 - tmp2) + A22 * (tmp1 - tmp2);
            double D1 = A12 * (tmp4 - tmp3) + A23 * (tmp2 - tmp1) + 2 * A13 * (tmp2 - tmp3);
            double D2 = A11 * (tmp2 - tmp3) + A12 * (tmp2 - tmp1);
            double D3 = A33 * (tmp2 - tmp3) + A23 * (tmp4 - tmp3);

            tmp1 = C1 * v1 [k] + D1 * v2 [k] + C2 * v3 [k];
            tmp2 = C2 * v1 [k] + D2 * v2 [k];
            tmp3 = C1 * v3 [k] + D3 * v2 [k];

            C1 = tmp1 * mn_norm + tmp2 * m_norm + tmp3 * n_norm;

            tmp1 = sinus * C1;

            if (k == 0) coord_vel_force_struct.forcex [I [l]] += tmp1;
            if (k == 1) coord_vel_force_struct.forcey [I [l]] += tmp1;
            if (k == 2) coord_vel_force_struct.forcez [I [l]] += tmp1;
        }
    }
}

////////////////////////////////////////////////////////////////////////
void polymer_constraints (double * box,
                          const std::vector<_Res_cg> &res,
                          const _topol_struct     &topol_struct,
                          _coord_vel_force_struct &coord_vel_force_struct,
                          _energy_struct          &energy_struct) {

    /////// RI - base
    for (int entry = 1; entry <= topol_struct.Nres_biomol; entry++) {
        int i = topol_struct.atom_key [entry][1];
        int j = topol_struct.atom_key [entry][2];
        if (i == 0 || j == 0) continue;

        int RES = topol_struct.res_simu2res_pdb[entry];
        std::string resname = res [entry - 1].name;

        if (resname == "A") {
            gen_bond_force (i, j, 10.0, 4.8515, coord_vel_force_struct, energy_struct.E_BOND);
            //bead_bond_force (i, j, 10.0, 4.8515, box, topol_struct, coord_vel_force_struct, energy_struct.E_BOND);

            // Base - RI - PH angle
            if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES) && \
                    entry < topol_struct.Nres_biomol) {
                int k = topol_struct.atom_key [entry + 1][0];
                gen_valence_force (j, i, k, 5.0, 1.9259, coord_vel_force_struct, energy_struct.E_ANGLE);
                //bead_valence_force (j, i, k, 5.0, 1.9259, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);
            }

            int k = topol_struct.atom_key [entry][0];
            if (k > 0) gen_valence_force (j, i, k, 5.0, 1.7029, coord_vel_force_struct, energy_struct.E_ANGLE);
                //bead_valence_force (j, i, k, 5.0, 1.7029, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);

        } else if (resname == "C") {
            gen_bond_force (i, j, 10.0, 4.2738, coord_vel_force_struct, energy_struct.E_BOND);
            //bead_bond_force (i, j, 10.0, 4.2738, box, topol_struct, coord_vel_force_struct, energy_struct.E_BOND);

            if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES) && \
                    entry < topol_struct.Nres_biomol) {
                int k = topol_struct.atom_key [entry + 1][0];
                gen_valence_force (j, i, k, 5.0, 1.9655, coord_vel_force_struct, energy_struct.E_ANGLE);
                //bead_valence_force (j, i, k, 5.0, 1.9655, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);
            }

            int k = topol_struct.atom_key [entry][0];
            if (k > 0) gen_valence_force (j, i, k, 5.0, 1.5803, coord_vel_force_struct, energy_struct.E_ANGLE);
                //bead_valence_force (j, i, k, 5.0, 1.5803, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);

        } else if (resname == "G") {
            gen_bond_force (i, j, 10.0, 4.9659, coord_vel_force_struct, energy_struct.E_BOND);
            //bead_bond_force (i, j, 10.0, 4.9659, box, topol_struct, coord_vel_force_struct, energy_struct.E_BOND);

            if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES) && \
                    entry < topol_struct.Nres_biomol) {
                int k = topol_struct.atom_key [entry + 1][0];
                gen_valence_force (j, i, k, 5.0, 1.9150, coord_vel_force_struct, energy_struct.E_ANGLE);
                //bead_valence_force (j, i, k, 5.0, 1.9150, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);
            }

            int k = topol_struct.atom_key [entry][0];
            if (k > 0) gen_valence_force (j, i, k, 5.0, 1.7690, coord_vel_force_struct, energy_struct.E_ANGLE);
                //bead_valence_force (j, i, k, 5.0, 1.7690, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);

        } else if (resname == "U") {
            gen_bond_force (i, j, 10.0, 4.2733, coord_vel_force_struct, energy_struct.E_BOND);
            //bead_bond_force (i, j, 10.0, 4.2733, box, topol_struct, coord_vel_force_struct, energy_struct.E_BOND);

            if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES) && \
                    entry < topol_struct.Nres_biomol) {
                int k = topol_struct.atom_key [entry + 1][0];
                gen_valence_force (j, i, k, 5.0, 1.9663, coord_vel_force_struct, energy_struct.E_ANGLE);
                //bead_valence_force (j, i, k, 5.0, 1.9663, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);
            }

            int k = topol_struct.atom_key [entry][0];
            if (k > 0) gen_valence_force (j, i, k, 5.0, 1.5735, coord_vel_force_struct, energy_struct.E_ANGLE);
                //bead_valence_force (j, i, k, 5.0, 1.5735, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);

        }
    }

    ////////// This goes downstream direction
    /// PH(n) - RI(n-1)

    for (int entry = 2; entry <= topol_struct.Nres_biomol; entry++) {
        int i = topol_struct.atom_key [entry][0];
        if (i == 0) continue;
        int prevRES = topol_struct.res_simu2res_pdb [entry - 1];
        if (std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), prevRES)) continue;
        int j = topol_struct.atom_key [entry - 1][1];
        if (j == 0) continue;

        gen_bond_force (i, j, 64.0, 3.8157, coord_vel_force_struct, energy_struct.E_BOND);
        //bead_bond_force (i, j, 64.0, 3.8157, box, topol_struct, coord_vel_force_struct, energy_struct.E_BOND);

        // PH(n) - RI(n-1) - PH(n-1) angle
        int k = topol_struct.atom_key [entry - 1][0];
        if (k > 0) gen_valence_force (i, j, k, 20.0, 1.4440, coord_vel_force_struct, energy_struct.E_ANGLE);
            //bead_valence_force (i, j, k, 20.0, 1.4440, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);

    }

    /// RI(n) - PH(n)

    for (int entry = 1; entry <= topol_struct.Nres_biomol; entry++) {
        int i = topol_struct.atom_key [entry][1];
        int j = topol_struct.atom_key [entry][0];
        if (i == 0 || j == 0) continue;

        gen_bond_force (i, j, 23.0, 4.6010, coord_vel_force_struct, energy_struct.E_BOND);
        //bead_bond_force (i, j, 23.0, 4.6010, box, topol_struct, coord_vel_force_struct, energy_struct.E_BOND);


        // RI(n) - PH(n) - RI(n-1) angle
        if (entry > 1) {
            int prevRES = topol_struct.res_simu2res_pdb [entry - 1];
            if (entry > 1 && not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), prevRES)) {
                int k = topol_struct.atom_key [entry - 1][1];
                if (k > 0) gen_valence_force (i, j, k, 20.0, 1.5256, coord_vel_force_struct, energy_struct.E_ANGLE);
                    //bead_valence_force (i, j, k, 20.0, 1.5256, box, topol_struct, coord_vel_force_struct, energy_struct.E_ANGLE);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////
void full_forces (const _topol_struct     &topol_struct,
                  _coord_vel_force_struct &coord_vel_force_struct,
                  unsigned int *myseed) {
    ////////// Construction of maxwell force, fill it with random numbers
//    for (int i = 1; i < 4501; i += 2) {
    for (size_t i = 1; i < coord_vel_force_struct.maxwell_force.size(); i += 2) {
        double r1;
#ifdef _OPENMP
        do r1 = (double) rand_r (myseed) / RAND_MAX;
#else
        do r1 = (double) rand () / RAND_MAX;
#endif
        while (r1 < 1e-9);

#ifdef _OPENMP
        double r2 = (double) rand_r (myseed) / RAND_MAX;
#else
        double r2 = (double) rand () / RAND_MAX;
#endif

        double tmp1 = sqrt (-2*log (r1));
        double tmp2 = 2.*PI * r2;

        coord_vel_force_struct.maxwell_force [i  ] = tmp1 * cos (tmp2);
        coord_vel_force_struct.maxwell_force [i+1] = tmp1 * sin (tmp2);
    }

    for (int i = 1; i <= topol_struct.Natm; i++) {
        int j = topol_struct.maxi_key [i];
        //int a_third = topol_struct.Natm / 3;
        //int k = i % a_third;

        coord_vel_force_struct.forcex [i] += coord_vel_force_struct.maxwell_force [i] * topol_struct.SIGMA_FORCE [j];
        coord_vel_force_struct.forcey [i] += coord_vel_force_struct.maxwell_force [i + topol_struct.Natm] * topol_struct.SIGMA_FORCE [j];
        coord_vel_force_struct.forcez [i] += coord_vel_force_struct.maxwell_force [i + 2*topol_struct.Natm] * topol_struct.SIGMA_FORCE [j];
    }
}

//////////////////////////////////////////////////////////////////////////////
// Apply constant force onto 2 pre-chosen beads (in the opposite directions)
//      (Force in pN, distance in Angstrom)
//////////////////////////////////////////////////////////////////////////////
void apply_external_force (int id1, int id2,
                           const _Dcoordinate &force, // direction of force applied onto bead id1
                                                      // the applied force onto bead id2 is in the opposite direction
                           _coord_vel_force_struct &coord_vel_force_struct) {
    coord_vel_force_struct.forcex [id1] += force.x;
    coord_vel_force_struct.forcey [id1] += force.y;
    coord_vel_force_struct.forcez [id1] += force.z;

    coord_vel_force_struct.forcex [id2] -= force.x;
    coord_vel_force_struct.forcey [id2] -= force.y;
    coord_vel_force_struct.forcez [id2] -= force.z;
}

#endif
