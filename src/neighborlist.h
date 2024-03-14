#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

//////// Separate particles in the simulation box to different cells,
///////          and count how many particles in each cell
void separate_in_cells (_simu_struct &simu_struct,
                        const _coord_vel_force_struct &coord_vel_force_struct,
                        int Natm) {
    //// Reset all
    if (simu_struct.Mc_total) {
        for (int m = 0; m <= simu_struct.Mc_total; m++)
            free (simu_struct.cell_content [m]);

        free (simu_struct.cell_content);
//        free (simu_struct.cell_mass);
//        simu_struct.cell_content.clear ();
        simu_struct.cell_mass.clear ();

        simu_struct.cell_content = NULL;
//        simu_struct.cell_mass = NULL;
    }

    // Divide simulation boxes into cells, separated by AT LEAST electrostatic cutoff distance
    simu_struct.Mc.x = (int) ceil (simu_struct.box[0] / simu_struct.R1_LIST);
    double lX = simu_struct.box[0] / simu_struct.Mc.x;

    simu_struct.Mc.y = (int) ceil (simu_struct.box[1] / simu_struct.R1_LIST);
    double lY = simu_struct.box[1] / simu_struct.Mc.y;

    simu_struct.Mc.z = (int) ceil (simu_struct.box[2] / simu_struct.R1_LIST);
    double lZ = simu_struct.box[2] / simu_struct.Mc.z;

    simu_struct.Mc_total = simu_struct.Mc.x * simu_struct.Mc.y * simu_struct.Mc.z;

    simu_struct.cell_content = (int **) calloc (simu_struct.Mc_total + 1, sizeof(int *));
//    for (int m = 0; m <= simu_struct.Mc_total; m++)
//        simu_struct.cell_content [m] = NULL;
//    simu_struct.cell_mass = (int *) calloc (simu_struct.Mc_total + 1, sizeof(int));

//    simu_struct.cell_content.resize (simu_struct.Mc_total + 1);
//    for (int i = 0; i < simu_struct.cell_content.size(); i++)
//        simu_struct.cell_content [i].push_back (0);       // not using 0 indices

    simu_struct.cell_mass.resize    (simu_struct.Mc_total + 1, 0);

    /////////////////////////////////////////////////////////////
    for (int i = 1; i <= Natm; i++) {
        double tmp = coord_vel_force_struct.coordx [i];
        if (tmp < 0)
            tmp += simu_struct.box[0];
        else if (tmp > simu_struct.box[0])
            tmp -= simu_struct.box[0];

        int i1 = (int) ceil (tmp / lX);

        tmp = coord_vel_force_struct.coordy [i];
        if (tmp < 0)
            tmp += simu_struct.box[1];
        else if (tmp > simu_struct.box[1])
            tmp -= simu_struct.box[1];

        int i2 = (int) ceil (tmp / lY);

        tmp = coord_vel_force_struct.coordz [i];
        if (tmp < 0)
            tmp += simu_struct.box[2];
        else if (tmp > simu_struct.box[2])
            tmp -= simu_struct.box[2];

        int i3 = (int) ceil (tmp / lZ);

        /////////////////////////////////////////////////////////////
        /// 1D Cell index
        int index = (i1 - 1) * simu_struct.Mc.y * simu_struct.Mc.z + (i2 - 1) * simu_struct.Mc.z + i3;

        simu_struct.cell_mass [index] ++;    // cell_mass stores how many particles in each cell

        simu_struct.cell_content [index] = (int *) realloc (simu_struct.cell_content [index], (simu_struct.cell_mass [index] + 1) * sizeof(int));
        simu_struct.cell_content [index] [simu_struct.cell_mass [index]] = i;    // cell_content stores which atoms (by the atom index) are in each cell
    }
}

//////////////////////////////////////////////////////////////////
void fill_neighbor_list (int i1, int i2,
                         _simu_struct                  &simu_struct,
                         const _topol_struct           &topol_struct,
                         const _coord_vel_force_struct &coord_vel_force_struct) {

    int m1 = topol_struct.maxi_key [i1];
    int m2 = topol_struct.maxi_key [i2];

    double vec [3];
    CC_vector (i1, i2, simu_struct.box, coord_vel_force_struct, vec);
    double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];

    double HB_cutoff = 9.0 + simu_struct.dr1;
    double HB2 = HB_cutoff*HB_cutoff;

    //// Ion - P case
    if ((m1 == 1 && m2 >= 7) || (m1 >= 7 && m2 == 1)) {
        if (std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), i1) || \
            std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), i2)) {
            if (r2 <= HB2) {
                simu_struct.HBneighborlist_mass += 2;
                simu_struct.HBneighborlist.push_back (i1);
                simu_struct.HBneighborlist.push_back (i2);
            }
        } else if (r2 <= simu_struct.R2_LIST) {
            simu_struct.DHneighborlist_mass += 2;
            simu_struct.DHneighborlist.push_back (i1);
            simu_struct.DHneighborlist.push_back (i2);
        }

    // Other cases
    } else {
        if (r2 <= HB2 && topol_struct.con_matrix [i1][i2] == 0) {
            simu_struct.HBneighborlist_mass += 2;
            simu_struct.HBneighborlist.push_back (i1);
            simu_struct.HBneighborlist.push_back (i2);
        }

        if (r2 <= simu_struct.R2_LIST && ((m1 == 1 && m2 == 1) || (m1 >= 7 && m2 >= 7)) && \
            not std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), i1) && \
            not std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), i2)) {
            simu_struct.DHneighborlist_mass += 2;
            simu_struct.DHneighborlist.push_back (i1);
            simu_struct.DHneighborlist.push_back (i2);
        }
    }

    // Electrostatic
    /*if (r2 <= simu_struct.R2_LIST && (m1 == 1 || m1 >= 7) && (m2 == 1 || m2 >= 7) && \
            not std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), i1) && \
            not std::binary_search (topol_struct.bead_zeroQ.begin(), topol_struct.bead_zeroQ.end(), i2)) {
        simu_struct.DHneighborlist_mass += 2;
        simu_struct.DHneighborlist.push_back (i1);
        simu_struct.DHneighborlist.push_back (i2);
    }

    // LJ, excluding ion-P
    // This list is built in order to monitor putative hydbonds formed during the course of simulation
    // 10A cutoff is considered safe
    int k;
    if (m1 < m2)
        k = (m1 - 1) * NA - (m1 - 1) * m1 / 2 + m2;
    else k = (m2 - 1) * NA - (m2 - 1) * m2 / 2 + m1;

    double HB_cutoff = 10. + simu_struct.dr1;

    if (r2 <= HB_cutoff*HB_cutoff && topol_struct.con_matrix [i1][i2] == 0 && \
            not ((m1 == 1 && m2 >= 7) || (m1 >= 7 && m2 == 1))) {
        simu_struct.HBneighborlist_mass += 2;
        simu_struct.HBneighborlist.push_back (i1);
        simu_struct.HBneighborlist.push_back (i2);
    }*/
}

//////////////////////////////////////////////////////////////////////////////////
/// Keep track of pairs of interaction (both LJ and electrostatic) within the CUTOFF
//////////////////////////////////////////////////////////////////////////////////
void populate_lists (_simu_struct                  &simu_struct,
                     const _topol_struct           &topol_struct,
                     const _coord_vel_force_struct &coord_vel_force_struct) {
    ///// index: Neighboring cells
    int index [13][3] = { { -1, -1, -1 }, { -1, -1, 0 }, { -1, -1, 1 },
                          { -1,  0, -1 }, { -1,  0, 0 }, { -1,  0, 1 },
                          { -1,  1, -1 }, { -1,  1, 0 }, { -1,  1, 1 },
                          {  0, -1, -1 }, {  0, -1, 0 }, {  0, -1, 1 }, { 0, 0, -1 } };

    separate_in_cells (simu_struct, coord_vel_force_struct, topol_struct.Natm);
//    simu_struct.list_content_SR.clear();
//    simu_struct.list_content_LR.clear();

//    simu_struct.list_mass_SR = 0;
//    simu_struct.list_mass_LR = 0;
    simu_struct.HBneighborlist_mass = 0;
    simu_struct.HBneighborlist.clear();

    simu_struct.DHneighborlist_mass = 0;
    simu_struct.DHneighborlist.clear();

//    simu_struct.list_content [0] = (int *) calloc (70000 + 1, sizeof(int));
//    simu_struct.list_content [1] = (int *) calloc (500000 + 1, sizeof(int));
//    simu_struct.list_content_SR.reserve (70000);       simu_struct.list_content_SR.push_back (0);
//    simu_struct.list_content_LR.reserve (500000);      simu_struct.list_content_LR.push_back (0);
    simu_struct.HBneighborlist.reserve (30000);        simu_struct.HBneighborlist.push_back (0);
    simu_struct.DHneighborlist.reserve (30000);        simu_struct.DHneighborlist.push_back (0);

    for (int m1 = 1; m1 <= simu_struct.Mc.x; m1++)
        for (int m2 = 1; m2 <= simu_struct.Mc.y; m2++)
            for (int m3 = 1; m3 <= simu_struct.Mc.z; m3++) {
                /// 1D index of the 'central' cell
                int J1 = (m1 - 1) * simu_struct.Mc.y * simu_struct.Mc.z + (m2 - 1) * simu_struct.Mc.z + m3;

                for (int k = 0; k < 13; k++) {
                    int n1 = m1 + index [k][0];

                    if (n1 > simu_struct.Mc.x) {
                        n1 -= simu_struct.Mc.x;
                        //if (n1 == m1 || n1 == m1 - 1) continue;
                    } else if (n1 < 1) {
                        n1 += simu_struct.Mc.x;
                        //if (n1 == m1 || n1 == m1 + 1) continue;
                    }

                    int n2 = m2 + index [k][1];
                    if (n2 > simu_struct.Mc.y) {
                        n2 -= simu_struct.Mc.y;
                        //if (n2 == m2 || n2 == m2 - 1) continue;
                    } else if (n2 < 1) {
                        n2 += simu_struct.Mc.y;
                        //if (n2 == m2 || n2 == m2 + 1) continue;
                    }

                    int n3 = m3 + index [k][2]; 
                    if (n3 > simu_struct.Mc.z) {
                        n3 -= simu_struct.Mc.z;
                        //if (n3 == m3 || n3 == m3 - 1) continue;
                    } else if (n3 < 1) {
                        n3 += simu_struct.Mc.z;
                        //if (n3 == m3 || n3 == m3 + 1) continue;
                    }

                    int J2 = (n1 - 1) * simu_struct.Mc.y * simu_struct.Mc.z + (n2 - 1) * simu_struct.Mc.z + n3;

                    for (int j1 = 1; j1 <= simu_struct.cell_mass [J1]; j1++) {
                        int i1 = simu_struct.cell_content [J1][j1];
//                        n1 = topol_struct.part_key [i1];

                        for (int j2 = 1; j2 <= simu_struct.cell_mass [J2]; j2++) {
                            int i2 = simu_struct.cell_content [J2][j2];
//                            n2 = topol_struct.part_key [i2];
                            fill_neighbor_list (i1, i2, simu_struct, topol_struct, coord_vel_force_struct);
                        }
                    }
                }

                ///////// Check interactions within this "central" cell
                for (int j1 = 1; j1 <= simu_struct.cell_mass [J1]; j1++) {
                    int i1 = simu_struct.cell_content [J1][j1];
                    // int n1 = topol_struct.part_key [i1];

                    for (int j2 = j1 + 1; j2 <= simu_struct.cell_mass [J1]; j2++) {
                        int i2 = simu_struct.cell_content [J1][j2];
                        // int n2 = topol_struct.part_key [i2];
                        fill_neighbor_list (i1, i2, simu_struct, topol_struct, coord_vel_force_struct);
                    }
                }
            }
}

//////////////////////////////////////////////////////////////
void check_shifts (_simu_struct                  &simu_struct,
                   const _topol_struct           &topol_struct,
                   _hydbond_struct               &hydbond_struct,
                   const _coord_vel_force_struct &coord_vel_force_struct) {

    double shift = 0;
    for (int i = 1; i <= topol_struct.Natm; i++) {
        double x = coord_vel_force_struct.coordx [i] - coord_vel_force_struct.coord_oldx [i];
        double y = coord_vel_force_struct.coordy [i] - coord_vel_force_struct.coord_oldy [i];
        double z = coord_vel_force_struct.coordz [i] - coord_vel_force_struct.coord_oldz [i];

        double r2 = x*x + y*y + z*z;
        if (r2 > shift) shift = r2;
    }

    shift = sqrt (shift);
    simu_struct.relist += (shift + shift);

//    if (simu_struct.relist > simu_struct.dr1 || simu_struct.relist_step == NSTEP_RELIST) {
    if (simu_struct.relist > simu_struct.dr1) {
        populate_lists (simu_struct, topol_struct, coord_vel_force_struct);
        if (simu_struct.fix_solute == 0 && hydbond_struct.eval) update_list_hydbond (simu_struct, topol_struct, hydbond_struct);
        simu_struct.relist = 0;
//        simu_struct.relist_step = 0;
    }
}

#endif
