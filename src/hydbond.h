#ifndef HYDBOND_H
#define HYDBOND_H

/// Update lists of possible hydbond interactions, because there is no need to do this every single step
void update_list_hydbond (const _simu_struct  &simu_struct,
                          const _topol_struct &topol_struct,
                          _hydbond_struct     &hydbond_struct) {
    for (int k = 1; k <= hydbond_struct.hb_N; k++)
        hydbond_struct.hb_status [k] = -1;

//    for (int k = 1; k <= simu_struct.list_mass_SR; k += 2) {
//        int i1 = simu_struct.list_content_SR [k];
    for (int k = 1; k <= simu_struct.HBneighborlist_mass; k += 2) {
        int i1 = simu_struct.HBneighborlist [k];
        if (topol_struct.part_key [i1] > 2) continue;

//        int i2 = simu_struct.list_content_SR [k + 1];
        int i2 = simu_struct.HBneighborlist [k+1];
        if (topol_struct.part_key [i2] > 2) continue;

        int J1;
        if (i1 < i2)
             J1 = (i1 - 1) * topol_struct.Natm_biomol - (i1 - 1) * i1 / 2 + i2;
        else J1 = (i2 - 1) * topol_struct.Natm_biomol - (i2 - 1) * i2 / 2 + i1;

        for (int j1 = 1; j1 <= hydbond_struct.HB_PAIR_N [J1]; j1++)
            hydbond_struct.hb_status [hydbond_struct.HB_PAIR [J1][j1]] = 0;
    }

    // Reset HB_EXCESS for every atomic atom
    // If the atom is not involved in the hydbond due to large distance, set HB_EXCESS less than 1 unit
    for (int k = 1; k <= hydbond_struct.HB_Natm; k++)
        hydbond_struct.HB_EXCESS [k] = hydbond_struct.HB_ATOM_N [k] - hydbond_struct.VALENCE [k];

    for (int k = 1; k <= hydbond_struct.hb_N; k++) {
        if (hydbond_struct.hb_status [k] == 0) continue;

        for (int i1 = 1; i1 <= hydbond_struct.ATOM_HB_N [k]; i1++)
            hydbond_struct.HB_EXCESS [hydbond_struct.ATOM_HB [k][i1]] -= 1;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////
void initialize_HB_atoms (const _Mol_cg       &mol_cg,
                          const _topol_struct &topol_struct,
                          _hydbond_struct     &hydbond_struct) {
//    hydbond_struct.HB_Nsite    = (int *)  calloc (topol_struct.Nres_biomol + 1, sizeof(int));
    hydbond_struct.HB_atom_key = (int **) calloc (topol_struct.Nres_biomol + 1, sizeof(int *));

    hydbond_struct.HB_Nsite.resize    (topol_struct.Nres_biomol + 1, 0);
//    hydbond_struct.HB_INDX.push_back (-1); //// Not using 0 indices
//    hydbond_struct.HB_JNDX.push_back (-1); ////      ''
    hydbond_struct.VALENCE.push_back (-1); ////      ''

//       0    1    2    3    4    5   6   7   8   9   10  11
//  A:  OP1  OP2  O2'  O3'  O4'  O5'  N1  N3  N6  N7  N9
//  C:  OP1  OP2  O2'  O3'  O4'  O5'  N1  O2  N3  N4
//  G:  OP1  OP2  O2'  O3'  O4'  O5'  N1  N2  N3  O6  N7  N9
//  U:  OP1  OP2  O2'  O3'  O4'  O5'  N1  O2  N3  O4

    for (int i = 1; i <= topol_struct.Nres_biomol; i++) {
        int k, l, m;  // m = -1 -> not count
        if (mol_cg.res[i-1].name == "A") {
            hydbond_struct.HB_Nsite [i] = 10;
            k = 2;   l = 8;   m = -1; // O2'   N6
        } else if (mol_cg.res[i-1].name == "C") {
            hydbond_struct.HB_Nsite [i] = 9;
            k = 2;   l = 7;   m = 9;  // O2'   O2   N4
        } else if (mol_cg.res[i-1].name == "G") {
            hydbond_struct.HB_Nsite [i] = 11;
            k = 2;   l = 7;   m = 9;  // O2'   N2   O6
        } else if (mol_cg.res[i-1].name == "U") {
            hydbond_struct.HB_Nsite [i] = 9;
            k = 2;   l = 7;   m = -1; // O2'   O2
        }
        hydbond_struct.HB_atom_key [i] = (int *) calloc (hydbond_struct.HB_Nsite [i] + 1, sizeof(int));

        for (int j = 0; j <= hydbond_struct.HB_Nsite [i]; j++) {
            hydbond_struct.HB_Natm ++;
            hydbond_struct.HB_atom_key [i][j] = hydbond_struct.HB_Natm;
            if (j == k || j == l || j == m)
                 hydbond_struct.VALENCE.push_back (2);
            else hydbond_struct.VALENCE.push_back (1);
        }
    }

//    hydbond_struct.HB_EXCESS = (int *)  calloc (hydbond_struct.HB_Natm + 1, sizeof(int));
//    hydbond_struct.HB_ATOM_N = (int *)  calloc (hydbond_struct.HB_Natm + 1, sizeof(int));
    hydbond_struct.HB_ATOM   = (int **) calloc (hydbond_struct.HB_Natm + 1, sizeof(int *));
//    hydbond_struct.EXCESS    = (int *)  calloc (hydbond_struct.HB_Natm + 1, sizeof(int));
//    hydbond_struct.HB_A      = (int *)  calloc (hydbond_struct.HB_Natm + 1, sizeof(int));

    hydbond_struct.HB_EXCESS.resize (hydbond_struct.HB_Natm + 1, 0);
    hydbond_struct.HB_ATOM_N.resize (hydbond_struct.HB_Natm + 1, 0);
//    hydbond_struct.HB_ATOM.resize   (hydbond_struct.HB_Natm + 1);
//    for (int i = 0; i < hydbond_struct.HB_ATOM.size(); i++)
//        hydbond_struct.HB_ATOM [i].push_back (0);

    hydbond_struct.EXCESS.resize    (hydbond_struct.HB_Natm + 1, 0);
    hydbond_struct.HB_A.resize      (hydbond_struct.HB_Natm + 1, 0);
}

//////////////////////////////////////////////////////////////////
int return_HB_JNDX (char * ATOM, const std::string &NUCLEOTIDE) {
    int j = -1;
    if      (strcmp (ATOM, "OP1") == 0) j = 0;
    else if (strcmp (ATOM, "OP2") == 0) j = 1;
    else if (strcmp (ATOM, "O2'") == 0) j = 2;
    else if (strcmp (ATOM, "O3'") == 0) j = 3;
    else if (strcmp (ATOM, "O4'") == 0) j = 4;
    else if (strcmp (ATOM, "O5'") == 0) j = 5;
    else if (strcmp (ATOM, "N1")  == 0) j = 6;
    else if (strcmp (ATOM, "O2")  == 0 || strcmp (ATOM, "N2") == 0) j = 7;
    else if (strcmp (ATOM, "N3")  == 0)
        if (NUCLEOTIDE == "A") j = 7;
        else j = 8;
    else if (strcmp (ATOM, "N6")  == 0) j = 8;
    else if (strcmp (ATOM, "O4")  == 0 || strcmp (ATOM, "O6") == 0 || strcmp (ATOM, "N4") == 0) j = 9;
    else if (strcmp (ATOM, "N7") == 0) {
        if      (NUCLEOTIDE == "A") j = 9;
        else if (NUCLEOTIDE == "G") j = 10;
    } else if (strcmp (ATOM, "N9") == 0) {
        if      (NUCLEOTIDE == "A") j = 10;
        else if (NUCLEOTIDE == "G") j = 11;
    }

    return j;
}

////////////////////////////////////////////////////////////////////////
/*int return_atom_key (bool phosphate_first, int resid,
                     const _topol_struct &topol_struct,
                     const std::string &group) {
    int atm_key;
    if (phosphate_first) {
        if      (group == "PH") atm_key = topol_struct.atom_key [2*resid - 1][0];
        else if (group == "RI") atm_key = topol_struct.atom_key [2*resid    ][0];
        else                    atm_key = topol_struct.atom_key [2*resid    ][1];
    } else {
        if      (group == "PH") atm_key = topol_struct.atom_key [2*resid - 2][0];
        else if (group == "RI") atm_key = topol_struct.atom_key [2*resid - 1][0];
        else                    atm_key = topol_struct.atom_key [2*resid - 1][1];
    }
    return atm_key;
}*/

///////////////////////////////////////////////////////////////////////////////////////////
void initialize_unprocessed_bonds (char * file_unprocessed_bonds,
                                   const _Mol_cg   &mol_cg,
                                   _topol_struct   &topol_struct,
                                   _hydbond_struct &hydbond_struct,
                                   const _coord_vel_force_struct &coord_vel_force_struct) {
    hydbond_struct.k_r     = 5.0;
    hydbond_struct.k_theta = 1.5;
    hydbond_struct.k_phi   = 0.15;

    FILE * f1;
    f1 = fopen (file_unprocessed_bonds, "r");
    fseek (f1, 0L, SEEK_END);
    long finish = ftell (f1);

    fseek (f1, 0L, SEEK_SET);
    int bp_N = 0;
    long position;
    while (position < finish) {
        bp_N ++;
        char A1 [10];

        do fscanf (f1, "%s", A1);
        while (strcmp (A1, "TER") != 0);
        position = ftell (f1);
    }

    ////////// Not using 0 indices
    hydbond_struct.hb_E.push_back   (0);
    hydbond_struct.hb_k10.push_back (0);
    hydbond_struct.hb_k11.push_back (0);
    hydbond_struct.hb_k12.push_back (0);
    hydbond_struct.hb_k20.push_back (0);
    hydbond_struct.hb_k21.push_back (0);
    hydbond_struct.hb_k22.push_back (0);

    hydbond_struct.hb_r.push_back      (0);
    hydbond_struct.hb_theta1.push_back (0);
    hydbond_struct.hb_theta2.push_back (0);
    hydbond_struct.hb_psi.push_back    (0);
    hydbond_struct.hb_psi1.push_back   (0);
    hydbond_struct.hb_psi2.push_back   (0);

    hydbond_struct.ATOM_HB_N.push_back (0);
    hydbond_struct.hb_dode.push_back   (-1);
    hydbond_struct.resid.push_back     (std::make_pair (-1, -1));
    hydbond_struct.ligand_bond.push_back (0);
//    hydbond_struct.hb_RES1.push_back   (0);
//    hydbond_struct.hb_RES2.push_back   (0);
    ///////////////////////////////////

    fseek (f1, 0L, SEEK_SET);
//    int status = 0;
    for (int i = 1; i < bp_N; i++) {
        int k_min = hydbond_struct.hb_N;
        char A0 [10], A01 [10], A02 [10];

        fscanf (f1, "%s", A0);
        fscanf (f1, "%s", A01);
        int RES1 = atoi (A01);
        int res1 = topol_struct.res_pdb2res_simu [RES1];
//        int RES1 = atoi (A01) - topol_struct.res_shift + 1;

        fscanf (f1, "%s", A02);
        int RES2 = atoi (A02);
        int res2 = topol_struct.res_pdb2res_simu [RES2];
//        int RES2 = atoi (A02) - topol_struct.res_shift + 1;

//        fscanf (f1, "%s", A1);
        char A1 [10] = " ", A2 [10] = " ";

        while (strcmp (A1, "TER") != 0) {
            fscanf (f1, "%s", A1);
            if (strcmp (A1, "TER") == 0) continue;
            fscanf (f1, "%s", A2);
            int k10, k11, k12, k20, k21, k22;
            bool ligand_bond = 0;

            // Pick relevant atoms for calculation of hbond geometry
            if (strcmp (A1, "OP1") == 0 || strcmp (A1, "OP2") == 0) {
//                k10 = return_atom_key (topol_struct.phosphate_first, RES1,     topol_struct, "RI");
//                k11 = return_atom_key (topol_struct.phosphate_first, RES1,     topol_struct, "PH");
//                k12 = return_atom_key (topol_struct.phosphate_first, RES1 + 1, topol_struct, "PH");
                k10 = topol_struct.atom_key [res1][1];
                k11 = topol_struct.atom_key [res1][0];
//                k10 = topol_struct.res2atmkey [RES1][1];
//                k11 = topol_struct.res2atmkey [RES1][0];
                if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES1) && \
//                    RES1 < topol_struct.res2atmkey.size() - 1)
                    res1 < topol_struct.Nres_biomol)
//                    k12 = topol_struct.res2atmkey [RES1+1][0];
                    k12 = topol_struct.atom_key [res1+1][0];
                else {
                    if (std::find (topol_struct.hbond_ignore.begin(), topol_struct.hbond_ignore.end(), RES1) == topol_struct.hbond_ignore.end())
                        topol_struct.hbond_ignore.push_back (RES1);
                    continue;
                }
            } else if (strcmp (A1, "O2'") == 0 || strcmp (A1, "O4'") == 0 || strcmp (A1, "O3'") == 0 || strcmp (A1, "O5'") == 0) {
//                k10 = return_atom_key (topol_struct.phosphate_first, RES1 + 1, topol_struct, "PH");
//                k11 = return_atom_key (topol_struct.phosphate_first, RES1,     topol_struct, "RI");
//                k12 = return_atom_key (topol_struct.phosphate_first, RES1 + 1, topol_struct, "RI");
                k11 = topol_struct.atom_key [res1][1];
//                k11 = topol_struct.res2atmkey [RES1][1];
                if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES1) && \
//                    RES1 < topol_struct.res2atmkey.size() - 1) {
                    res1 < topol_struct.Nres_biomol) {
//                    k10 = topol_struct.res2atmkey [RES1+1][0];
//                    k12 = topol_struct.res2atmkey [RES1+1][1];
                    k10 = topol_struct.atom_key [res1+1][0];
                    k12 = topol_struct.atom_key [res1+1][1];
                    if (k10 == 0 || k12 == 0) {
                        if (std::find (topol_struct.hbond_ignore.begin(), topol_struct.hbond_ignore.end(), RES1) == topol_struct.hbond_ignore.end())
                            topol_struct.hbond_ignore.push_back (RES1);
                        continue;
                    }
                } else {
                    if (std::find (topol_struct.hbond_ignore.begin(), topol_struct.hbond_ignore.end(), RES1) == topol_struct.hbond_ignore.end())
                        topol_struct.hbond_ignore.push_back (RES1);
                    continue;
                }
            } else {
//                k10 = return_atom_key (topol_struct.phosphate_first, RES1,     topol_struct, "RI");
//                k11 = return_atom_key (topol_struct.phosphate_first, RES1,     topol_struct, "Base");
//                k12 = return_atom_key (topol_struct.phosphate_first, RES1 + 1, topol_struct, "PH");
                k10 = topol_struct.atom_key [res1][1];
                k11 = topol_struct.atom_key [res1][2];
//                k10 = topol_struct.res2atmkey [RES1][1];
//                k11 = topol_struct.res2atmkey [RES1][2];
                if (RES1 == topol_struct.ligand_index) {
                    k12 = 0;
                    ligand_bond = 1;
                } else if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES1) && \
//                    RES1 < topol_struct.res2atmkey.size() - 1)
                    res1 < topol_struct.Nres_biomol)
//                    k12 = topol_struct.res2atmkey [RES1+1][0];
                    k12 = topol_struct.atom_key [res1+1][0];
                else {
                    if (std::find (topol_struct.hbond_ignore.begin(), topol_struct.hbond_ignore.end(), RES1) == topol_struct.hbond_ignore.end())
                        topol_struct.hbond_ignore.push_back (RES1);
                    continue;
                }
            }

            if (strcmp (A2, "OP1") == 0 || strcmp (A2, "OP2") == 0) { 
//                k20 = return_atom_key (topol_struct.phosphate_first, RES2,     topol_struct, "RI");
//                k21 = return_atom_key (topol_struct.phosphate_first, RES2,     topol_struct, "PH");
//                k22 = return_atom_key (topol_struct.phosphate_first, RES2 + 1, topol_struct, "PH");
                k20 = topol_struct.atom_key [res2][1];
                k21 = topol_struct.atom_key [res2][0];
//                k20 = topol_struct.res2atmkey [RES2][1];
//                k21 = topol_struct.res2atmkey [RES2][0];
                if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES2) && \
//                    RES2 < topol_struct.res2atmkey.size() - 1)
                    res2 < topol_struct.Nres_biomol)
//                    k22 = topol_struct.res2atmkey [RES2+1][0];
                    k22 = topol_struct.atom_key [res2+1][0];
                else {
                    if (std::find (topol_struct.hbond_ignore.begin(), topol_struct.hbond_ignore.end(), RES2) == topol_struct.hbond_ignore.end())
                        topol_struct.hbond_ignore.push_back (RES2);
                    continue;
                }

            } else if (strcmp (A2, "O2'") == 0 || strcmp (A2, "O4'") == 0 || strcmp (A2, "O3'") == 0 || strcmp (A2, "O5'") == 0) {
//                k20 = return_atom_key (topol_struct.phosphate_first, RES2 + 1, topol_struct, "PH");
//                k21 = return_atom_key (topol_struct.phosphate_first, RES2,     topol_struct, "RI");
//                k22 = return_atom_key (topol_struct.phosphate_first, RES2 + 1, topol_struct, "RI");
                k21 = topol_struct.atom_key [res2][1];
//                k21 = topol_struct.res2atmkey [RES2][1];
                if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES2) && \
//                    RES2 < topol_struct.res2atmkey.size() - 1) {
                    res2 < topol_struct.Nres_biomol) {
                    k20 = topol_struct.atom_key [res2+1][0];
                    k22 = topol_struct.atom_key [res2+1][1];
//                    k20 = topol_struct.res2atmkey [RES2+1][0];
//                    k22 = topol_struct.res2atmkey [RES2+1][1];
                    if (k20 == 0 || k22 == 0) {
                        if (std::find (topol_struct.hbond_ignore.begin(), topol_struct.hbond_ignore.end(), RES2) == topol_struct.hbond_ignore.end())
                            topol_struct.hbond_ignore.push_back (RES2);
                        continue;
                    }
                } else {
                    if (std::find (topol_struct.hbond_ignore.begin(), topol_struct.hbond_ignore.end(), RES2) == topol_struct.hbond_ignore.end())
                        topol_struct.hbond_ignore.push_back (RES2);
                    continue;
                }

            } else {
//                k20 = return_atom_key (topol_struct.phosphate_first, RES2,     topol_struct, "RI");
//                k21 = return_atom_key (topol_struct.phosphate_first, RES2,     topol_struct, "Base");
//                k22 = return_atom_key (topol_struct.phosphate_first, RES2 + 1, topol_struct, "PH");
                k20 = topol_struct.atom_key [res2][1];
                k21 = topol_struct.atom_key [res2][2];
//                k20 = topol_struct.res2atmkey [RES2][1];
//                k21 = topol_struct.res2atmkey [RES2][2];
                if (RES2 == topol_struct.ligand_index) {
                    k22 = 0;
                    ligand_bond = 1;
                } else if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES2) && \
//                    RES2 < topol_struct.res2atmkey.size() - 1)
                    res2 < topol_struct.Nres_biomol)
//                    k22 = topol_struct.res2atmkey [RES2+1][0];
                    k22 = topol_struct.atom_key [res2+1][0];
                else {
                    if (std::find (topol_struct.hbond_ignore.begin(), topol_struct.hbond_ignore.end(), RES2) == topol_struct.hbond_ignore.end())
                        topol_struct.hbond_ignore.push_back (RES2);
                    continue;
                }
            }

            // Check if this bond had been initialized
            bool status = 0;
            int k;
            for (k = k_min + 1; k <= hydbond_struct.hb_N; k++)
                if (k10 == hydbond_struct.hb_k10 [k] && \
                    k11 == hydbond_struct.hb_k11 [k] && \
                    k12 == hydbond_struct.hb_k12 [k] && \
                    k20 == hydbond_struct.hb_k20 [k] && \
                    k21 == hydbond_struct.hb_k21 [k] && \
                    k22 == hydbond_struct.hb_k22 [k]) {
                    status = 1;
                    break;
                }

            if (status) {
                hydbond_struct.hb_E [k] -= hydbond_struct.hs_D;// * 0.73753;
                status = 0;
           /* } else if (strcmp (A0, "CAN51") == 0 && strcmp (A2, "O2'") == 0 \
                    || strcmp (A0, "CAN52") == 0 && strcmp (A1, "O2'") == 0) {
                k = hydbond_struct.hb_N;*/
                //hb_E [k] -= hs_D * 0.73753;
            } else {
                hydbond_struct.hb_N++;
                hydbond_struct.hb_k10.push_back (k10);
                hydbond_struct.hb_k11.push_back (k11);
                hydbond_struct.hb_k12.push_back (k12);
                hydbond_struct.hb_k20.push_back (k20);
                hydbond_struct.hb_k21.push_back (k21);
                hydbond_struct.hb_k22.push_back (k22);
                hydbond_struct.ligand_bond.push_back (ligand_bond);

                hydbond_struct.hb_E.push_back      (- hydbond_struct.hs_D);// * 0.73753);

                // Check those HARD-CODE
                if (strcmp (A0, "CAN11") == 0 || strcmp (A0, "NATA-U") == 0 || strcmp (A0, "A-U") == 0) {          // A - U
                    hydbond_struct.hb_r     .push_back (5.8815);
                    hydbond_struct.hb_theta1.push_back (2.7283);
                    hydbond_struct.hb_theta2.push_back (2.5117);
                    hydbond_struct.hb_psi   .push_back (1.2559);
                    hydbond_struct.hb_psi1  .push_back (0.9545);
                    hydbond_struct.hb_psi2  .push_back (1.1747);

                } else if (strcmp (A0, "CAN12") == 0 || strcmp (A0, "NATU-A") == 0 || strcmp (A0, "U-A") == 0) {   // U - A
                    hydbond_struct.hb_r     .push_back (5.8815);
                    hydbond_struct.hb_theta1.push_back (2.5117);
                    hydbond_struct.hb_theta2.push_back (2.7283);
                    hydbond_struct.hb_psi   .push_back (1.2559);
                    hydbond_struct.hb_psi1  .push_back (1.1747);
                    hydbond_struct.hb_psi2  .push_back (0.9545);

                } else if (strcmp (A0, "CAN21") == 0 || strcmp (A0, "NATC-G") == 0 || strcmp (A0, "C-G") == 0) {   // C - G
                    hydbond_struct.hb_r     .push_back (5.6550);
                    hydbond_struct.hb_theta1.push_back (2.4837);
                    hydbond_struct.hb_theta2.push_back (2.8230);
                    hydbond_struct.hb_psi   .push_back (1.3902);
                    hydbond_struct.hb_psi1  .push_back (1.2174);
                    hydbond_struct.hb_psi2  .push_back (0.7619);

                } else if (strcmp (A0, "CAN22") == 0 || strcmp (A0, "NATG-C") == 0 || strcmp (A0, "G-C") == 0) {   // G - C
                    hydbond_struct.hb_r     .push_back (5.6550);
                    hydbond_struct.hb_theta1.push_back (2.8230);
                    hydbond_struct.hb_theta2.push_back (2.4837);
                    hydbond_struct.hb_psi   .push_back (1.3902);
                    hydbond_struct.hb_psi1  .push_back (0.7619);
                    hydbond_struct.hb_psi2  .push_back (1.2174);

                } else if (strcmp (A0, "CAN31") == 0 || strcmp (A0, "NATG-U") == 0 || strcmp (A0, "G-U") == 0) {   // G - U
                    hydbond_struct.hb_r     .push_back (5.9234);
                    hydbond_struct.hb_theta1.push_back (2.9051);
                    hydbond_struct.hb_theta2.push_back (2.0544);
                    hydbond_struct.hb_psi   .push_back (2.1469);
                    hydbond_struct.hb_psi1  .push_back (-0.5366);
                    hydbond_struct.hb_psi2  .push_back (1.4288);

                } else if (strcmp (A0, "CAN32") == 0 || strcmp (A0, "NATU-G") == 0 || strcmp (A0, "U-G") == 0) {   // U - G
                    hydbond_struct.hb_r     .push_back (5.9234);
                    hydbond_struct.hb_theta1.push_back (2.0544);
                    hydbond_struct.hb_theta2.push_back (2.9051);
                    hydbond_struct.hb_psi   .push_back (2.1469);
                    hydbond_struct.hb_psi1  .push_back (1.4288);
                    hydbond_struct.hb_psi2  .push_back (-0.5366);

                //////// What are CAN41, CAN42, CAN51 and CAN52?
                } else if (strcmp (A0, "CAN41") == 0) {
                    hydbond_struct.hb_r     .push_back (5.3737);
                    hydbond_struct.hb_theta1.push_back (2.8295);
                    hydbond_struct.hb_theta2.push_back (2.0314);
                    hydbond_struct.hb_psi   .push_back (0.6063);
                    hydbond_struct.hb_psi1  .push_back (1.0893);
                    hydbond_struct.hb_psi2  .push_back (1.4835);

                } else if (strcmp (A0, "CAN42") == 0) {
                    hydbond_struct.hb_r     .push_back (5.3737);
                    hydbond_struct.hb_theta1.push_back (2.0314);
                    hydbond_struct.hb_theta2.push_back (2.8295);
                    hydbond_struct.hb_psi   .push_back (0.6063);
                    hydbond_struct.hb_psi1  .push_back (1.4835);
                    hydbond_struct.hb_psi2  .push_back (1.0893);

                } else if (strcmp (A0, "CAN51") == 0) {
                    hydbond_struct.hb_r     .push_back (6.4922);
                    hydbond_struct.hb_theta1.push_back (2.0579);
                    hydbond_struct.hb_theta2.push_back (1.3569);
                    hydbond_struct.hb_psi   .push_back (3.1147);
                    hydbond_struct.hb_psi1  .push_back (-1.4796);
                    hydbond_struct.hb_psi2  .push_back (1.5278);

                } else if (strcmp (A0, "CAN52") == 0) {
                    hydbond_struct.hb_r     .push_back (6.4922);
                    hydbond_struct.hb_theta1.push_back (1.3569);
                    hydbond_struct.hb_theta2.push_back (2.0579);
                    hydbond_struct.hb_psi   .push_back (3.1147);
                    hydbond_struct.hb_psi1  .push_back (1.5278);
                    hydbond_struct.hb_psi2  .push_back (-1.4796);

                } else {
                    /// Tertiary hbond, push values found in crystal structure
                    hydbond_struct.hb_r     .push_back (CC_distance    (k11, k21, coord_vel_force_struct));
                    if (not hydbond_struct.ligand_bond.back()) {
                        hydbond_struct.hb_theta1.push_back (valence_angle  (k10, k11, k21, coord_vel_force_struct));
                        hydbond_struct.hb_theta2.push_back (valence_angle  (k20, k21, k11, coord_vel_force_struct));
                        hydbond_struct.hb_psi   .push_back (dihedral_angle (k10, k11, k21, k20, coord_vel_force_struct));
                        hydbond_struct.hb_psi1  .push_back (dihedral_angle (k21, k11, k10, k12, coord_vel_force_struct));
                        hydbond_struct.hb_psi2  .push_back (dihedral_angle (k11, k21, k20, k22, coord_vel_force_struct));
                    } else {
                        hydbond_struct.hb_theta1.push_back (0);
                        hydbond_struct.hb_theta2.push_back (0);
                        hydbond_struct.hb_psi   .push_back (0);
                        hydbond_struct.hb_psi1  .push_back (0);
                        hydbond_struct.hb_psi2  .push_back (0);
                    }
                }

                /////////////////////////////////////////////////////////////////
//                hydbond_struct.ATOM_HB_N = (int *) realloc (hydbond_struct.ATOM_HB_N, (hydbond_struct.hb_N + 1) * sizeof(int));
//                hydbond_struct.ATOM_HB_N [hydbond_struct.hb_N] = 0;

                hydbond_struct.ATOM_HB_N.push_back (0);  ////// initialize ATOM_HB_N [current_index]
                hydbond_struct.ATOM_HB = (int **) realloc (hydbond_struct.ATOM_HB, (hydbond_struct.hb_N + 1) * sizeof(int *));
                hydbond_struct.ATOM_HB [hydbond_struct.hb_N] = NULL;
//                hydbond_struct.ATOM_HB.resize (hydbond_struct.hb_N + 1);
//                hydbond_struct.ATOM_HB [hydbond_struct.hb_N].push_back (0);

//                hydbond_struct.hb_dode = (int *) realloc (hydbond_struct.hb_dode, (hydbond_struct.hb_N + 1) * sizeof(int));

//                if (strcmp (A0, "NON00") == 0 || strcmp (A0, "APLAT") == 0)
//                     hydbond_struct.hb_dode [hydbond_struct.hb_N] = 1;
//                else hydbond_struct.hb_dode [hydbond_struct.hb_N] = 0;

                if (strcmp (A0, "NON00") == 0 || strcmp (A0, "APLAT") == 0 || strstr (A0, "NAT") != NULL)
                     hydbond_struct.hb_dode.push_back (1);       // native     bond
                else hydbond_struct.hb_dode.push_back (0);       // non-native hbond

                hydbond_struct.hb_code = (char **) realloc (hydbond_struct.hb_code, (hydbond_struct.hb_N + 1) * sizeof(char *));
                hydbond_struct.hb_code [hydbond_struct.hb_N] = (char *) calloc (10, sizeof(char));

                strcpy (hydbond_struct.hb_code [hydbond_struct.hb_N], A0);

                hydbond_struct.resid.push_back (std::make_pair (res1, res2));

/*                hydbond_struct.hb_RES1 = (int *) realloc (hydbond_struct.hb_RES1, (hydbond_struct.hb_N + 1) * sizeof(int)); 
                hydbond_struct.hb_RES1 [hydbond_struct.hb_N] = atoi (A01);
                hydbond_struct.hb_RES2 = (int *) realloc (hydbond_struct.hb_RES2, (hydbond_struct.hb_N + 1) * sizeof(int));
                hydbond_struct.hb_RES2 [hydbond_struct.hb_N] = atoi (A02);*/

//                hydbond_struct.hb_RES1.push_back (atoi (A01));
//                hydbond_struct.hb_RES2.push_back (atoi (A02));
                k = hydbond_struct.hb_N;
            }

            // Keep track of which atoms (in the atomic structure) participating in hbond
            // HB_atom_key [residue][HB_JNDX] = HB_Natm
            //      |                             |
            //      v                             v
            // index of atom in hydbond      total No of atoms that **CAN** hbond (in the atomistic structure)
            //                 (see initialize_HB_atoms and HB_JNDX)
            
            // ATOM_HB_N [hb_N] : No of atoms involving in a single coarse_grained hydbond
            //             |
            //             v
            //         total No of considered coarse-grained hbonds

            // ATOM_HB [hb_N] : stores which atoms involving in hbond, base on HB_atom_key
            // HB_ATOM_N [HB_Natm] : how many hydbonds an (atomic) atom is involving
            // HB_ATOM [HB_Natm] : stores indices (1..hb_N) of coarse_grained hydbond an atom is involving
            // HB_EXCESS [HB_Natm] = HB_ATOM_N [HB_Natm] - VALENCE [HB_Natm]
            //         stores how many hydbonds need to be got rid of

            int hb_jndx = return_HB_JNDX (A1, mol_cg.res[res1 - 1].name);
            int hb_key1 = hydbond_struct.HB_atom_key [res1][hb_jndx];

            hb_jndx = return_HB_JNDX (A2, mol_cg.res[res2 - 1].name);
            int hb_key2 = hydbond_struct.HB_atom_key [res2][hb_jndx];

            hydbond_struct.ATOM_HB_N [k] += 2;
//            hydbond_struct.ATOM_HB = (int **) realloc (hydbond_struct.ATOM_HB, (k + 1) * sizeof(int *));
//            hydbond_struct.ATOM_HB [k] = NULL;
            hydbond_struct.ATOM_HB [k] = (int *) realloc (hydbond_struct.ATOM_HB [k], (hydbond_struct.ATOM_HB_N [k] + 1) * sizeof(int));
            hydbond_struct.ATOM_HB [k][hydbond_struct.ATOM_HB_N [k] - 1] = hb_key1;
            hydbond_struct.ATOM_HB [k][hydbond_struct.ATOM_HB_N [k]]     = hb_key2;
//            hydbond_struct.ATOM_HB [k].push_back (k12);
//            hydbond_struct.ATOM_HB [k].push_back (k22);

            hydbond_struct.HB_ATOM_N [hb_key1] ++;
            hydbond_struct.HB_ATOM [hb_key1] = (int *) realloc (hydbond_struct.HB_ATOM [hb_key1], (hydbond_struct.HB_ATOM_N [hb_key1] + 1) * sizeof(int));
            hydbond_struct.HB_ATOM [hb_key1][hydbond_struct.HB_ATOM_N [hb_key1]] = k;
//            hydbond_struct.HB_ATOM [k12].push_back (k);

            hydbond_struct.HB_ATOM_N [hb_key2] ++;
            hydbond_struct.HB_ATOM [hb_key2] = (int *) realloc (hydbond_struct.HB_ATOM [hb_key2], (hydbond_struct.HB_ATOM_N [hb_key2] + 1) * sizeof(int));
            hydbond_struct.HB_ATOM [hb_key2][hydbond_struct.HB_ATOM_N [hb_key2]] = k;
//            hydbond_struct.HB_ATOM [k22].push_back (k);
//            fscanf (f1, "%s", A1);
        } 
    }
    fclose (f1);

    if (topol_struct.hbond_ignore.size() > 0) {
        printf ("Ignoring Hbonds from those nucleotides: ");
        for (size_t i = 0; i < topol_struct.hbond_ignore.size(); i++)
            printf ("%d  ", topol_struct.hbond_ignore [i]);
        printf ("\n");
    }

    hydbond_struct.hb_psi0.resize   (hydbond_struct.hb_N + 1, 0);
    hydbond_struct.hb_psi10.resize  (hydbond_struct.hb_N + 1, 0);
    hydbond_struct.hb_psi20.resize  (hydbond_struct.hb_N + 1, 0);
    hydbond_struct.hb_energy.resize (hydbond_struct.hb_N + 1, 0);
    hydbond_struct.hb_f.resize      (hydbond_struct.hb_N + 1, 0);
    hydbond_struct.hb_status.resize (hydbond_struct.hb_N + 1, 0);
    hydbond_struct.hb_K.resize      (hydbond_struct.hb_N + 1, 0);

    topol_struct.Npair_biomol = topol_struct.Natm_biomol * (topol_struct.Natm_biomol + 1) / 2;
//    hydbond_struct.HB_PAIR_N = (int *)  calloc (topol_struct.Npair_biomol + 1, sizeof(int));
    hydbond_struct.HB_PAIR   = (int **) calloc (topol_struct.Npair_biomol + 1, sizeof(int *));
    hydbond_struct.HB_PAIR_N.resize (topol_struct.Npair_biomol + 1, 0);
//    hydbond_struct.HB_PAIR.resize   (topol_struct.Npair_biomol + 1);

    // HB_PAIR_N stores the number of hydbond a particular pair of beads participating
    // HB_PAIR points to the hydbond index in which two particular beads participate
    // hb_status:   -1 : distance between two beads > short range interaction cutoff -> hydbond ignored until the neighbor list updated
    //               1 : distance between two beads is not within 2A from the "standard" structure -> hydbond will be ignored
    //               0 : distance between two beads is within 2A from the "standard" structure -> hydbond will be calculated

    for (int i = 1; i <= hydbond_struct.hb_N; i++) {
        hydbond_struct.hb_status [i] = -1;

        int k11 = hydbond_struct.hb_k11 [i];
        int k21 = hydbond_struct.hb_k21 [i];

        int k;
        if (k11 < k21)
             k = (k11 - 1) * topol_struct.Natm_biomol - (k11 - 1) * k11 / 2 + k21;
        else k = (k21 - 1) * topol_struct.Natm_biomol - (k21 - 1) * k21 / 2 + k11;

        hydbond_struct.HB_PAIR_N [k] ++;
//        hydbond_struct.HB_PAIR [k].push_back (i);
        hydbond_struct.HB_PAIR [k] = (int *) realloc (hydbond_struct.HB_PAIR [k], (hydbond_struct.HB_PAIR_N [k] + 1) * sizeof(int));
        hydbond_struct.HB_PAIR [k][hydbond_struct.HB_PAIR_N [k]] = i;
    }
}

///////////////////////////////////////////////////////
void compute_hbond (_coord_vel_force_struct &coord_vel_force_struct,
                    _hydbond_struct         &hydbond_struct) {
    
    for (int i = 1; i <= hydbond_struct.HB_Natm; i++)
        hydbond_struct.EXCESS [i] = hydbond_struct.HB_EXCESS [i];

    for (int k = 1; k <= hydbond_struct.hb_N; k++) {
        if (hydbond_struct.hb_status [k] == -1) continue;

        // Ignore bonds that have large distance deviation
        double r = CC_distance (hydbond_struct.hb_k11 [k], hydbond_struct.hb_k21 [k], coord_vel_force_struct);
        //if ((fabs (r - hydbond_struct.hb_r [k]) >  2.0 && hydbond_struct.hb_dode [k] == 0) || \
            (fabs (r - hydbond_struct.hb_r [k]) > 10.0 && hydbond_struct.hb_dode [k] == 1)) {
        if (fabs (r - hydbond_struct.hb_r [k]) > 2.0) {
            hydbond_struct.hb_status [k] = 1;
            for (int j = 1; j <= hydbond_struct.ATOM_HB_N [k]; j++)
                hydbond_struct.EXCESS [hydbond_struct.ATOM_HB [k][j]] -= 1;
            continue;

        // Native hbond is more important than non-native ones
        // Once a native bond is formed, do not consider non-native ones
        }/* else if (hydbond_struct.hb_dode [k]) {
            int res1 = hydbond_struct.resid [k].first;
            int res2 = hydbond_struct.resid [k].second;
            for (int l = 1; l <= hydbond_struct.hb_N; l++) {
                if (hydbond_struct.hb_status [l] == -1 || hydbond_struct.hb_status [l] == 1 || hydbond_struct.hb_dode [l]) continue;
                if (hydbond_struct.resid [l].first == res1 || hydbond_struct.resid [l].second == res1 || \
                    hydbond_struct.resid [l].first == res2 || hydbond_struct.resid [l].second == res2) {
                    hydbond_struct.hb_status [l] = 1;
                    for (int j = 1; j <= hydbond_struct.ATOM_HB_N [l]; j++)
                        hydbond_struct.EXCESS [hydbond_struct.ATOM_HB [l][j]] -= 1;
                }
            }
        }*/

        double r2 = hydbond_struct.k_r * (r - hydbond_struct.hb_r [k]) * (r - hydbond_struct.hb_r [k]);

        if (not hydbond_struct.ligand_bond [k]) {
            double theta1 = valence_angle (hydbond_struct.hb_k10 [k], hydbond_struct.hb_k11 [k], hydbond_struct.hb_k21 [k], coord_vel_force_struct);
            double theta2 = valence_angle (hydbond_struct.hb_k20 [k], hydbond_struct.hb_k21 [k], hydbond_struct.hb_k11 [k], coord_vel_force_struct);

            double psi  = dihedral_angle (hydbond_struct.hb_k10 [k], hydbond_struct.hb_k11 [k], hydbond_struct.hb_k21 [k], hydbond_struct.hb_k20 [k], coord_vel_force_struct);
            double psi1 = dihedral_angle (hydbond_struct.hb_k21 [k], hydbond_struct.hb_k11 [k], hydbond_struct.hb_k10 [k], hydbond_struct.hb_k12 [k], coord_vel_force_struct);
            double psi2 = dihedral_angle (hydbond_struct.hb_k11 [k], hydbond_struct.hb_k21 [k], hydbond_struct.hb_k20 [k], hydbond_struct.hb_k22 [k], coord_vel_force_struct);

            ///////// Should reuse code from stacking_interactions????
            if (fabs (psi - hydbond_struct.hb_psi [k]) <= PI)
                hydbond_struct.hb_psi0 [k] = hydbond_struct.hb_psi [k];
            else if (psi - hydbond_struct.hb_psi [k] > PI)
                hydbond_struct.hb_psi0 [k] = hydbond_struct.hb_psi [k] + 2*PI;
            else hydbond_struct.hb_psi0 [k] = hydbond_struct.hb_psi [k] - 2*PI;

            if (fabs (psi1 - hydbond_struct.hb_psi1 [k]) <= PI)
                hydbond_struct.hb_psi10 [k] = hydbond_struct.hb_psi1 [k];
            else if (psi1 - hydbond_struct.hb_psi1 [k] > PI)
                hydbond_struct.hb_psi10 [k] = hydbond_struct.hb_psi1 [k] + 2*PI;
            else hydbond_struct.hb_psi10 [k] = hydbond_struct.hb_psi1 [k] - 2*PI;

            if (fabs (psi2 - hydbond_struct.hb_psi2 [k]) <= PI)
                hydbond_struct.hb_psi20 [k] = hydbond_struct.hb_psi2 [k];
            else if (psi2 - hydbond_struct.hb_psi2 [k] > PI)
                hydbond_struct.hb_psi20 [k] = hydbond_struct.hb_psi2 [k] + 2*PI;
            else hydbond_struct.hb_psi20 [k] = hydbond_struct.hb_psi2 [k] - 2*PI;

            /// Check HARD-CODE
            r2 += hydbond_struct.k_theta * (theta1 - hydbond_struct.hb_theta1 [k]) * (theta1 - hydbond_struct.hb_theta1 [k]);
            r2 += hydbond_struct.k_theta * (theta2 - hydbond_struct.hb_theta2 [k]) * (theta2 - hydbond_struct.hb_theta2 [k]);
            r2 += hydbond_struct.k_phi * (psi  - hydbond_struct.hb_psi0  [k]) * (psi  - hydbond_struct.hb_psi0  [k]);
            r2 += hydbond_struct.k_phi * (psi1 - hydbond_struct.hb_psi10 [k]) * (psi1 - hydbond_struct.hb_psi10 [k]);
            r2 += hydbond_struct.k_phi * (psi2 - hydbond_struct.hb_psi20 [k]) * (psi2 - hydbond_struct.hb_psi20 [k]);
        }

        /*if (hydbond_struct.hb_dode [k] == 0) {
            hydbond_struct.hb_energy [k] = hydbond_struct.hb_E [k] * exp (-r2);
            hydbond_struct.hb_f [k] = - hydbond_struct.hb_energy [k];
        } else {
            hydbond_struct.hb_energy [k] = hydbond_struct.hb_E [k] / (1. + r2);
            hydbond_struct.hb_f [k] = - hydbond_struct.hb_energy [k] / (1. + r2);
        }*/
        hydbond_struct.hb_energy [k] = hydbond_struct.hb_E [k] * exp(-r2);
        hydbond_struct.hb_f [k] = - hydbond_struct.hb_energy [k];

        /// Reduce non-native interaction weight
        //if (hydbond_struct.hb_dode [k] == 0)
        //    hydbond_struct.hb_energy [k] *= hydbond_struct.weight_nonnative;
    }
}

////////////////////////////////////////////////////////
void remove_excess_hbond (_hydbond_struct &hydbond_struct,
                          const double &T,
                          unsigned int *myseed) {
    // Removing excess Hbond, nA stores how many atoms that have excess Hbond
    //                        HB_A stores those atom locations
    //                        nK stores how many atomic Hbonds involving in this Monte-Carlo process
    //                        hb_K stores which coarse-grained beads involving
    while (true) {
        int nA = 0, nK = 0;
        for (int i = 1; i <= hydbond_struct.HB_Natm; i++)
            if (hydbond_struct.EXCESS [i] > 0) {
                nA += 1;
                hydbond_struct.HB_A [nA] = i;
            }

        if (nA == 0) break;

        // Randomly pick an atom
#ifdef _OPENMP
        double r = (double) rand_r (myseed) / RAND_MAX;
#else
        double r = (double) rand () / RAND_MAX;
#endif
        int tmp = (int) ceil (r * nA);
        if      (tmp == nA + 1) tmp = nA;
        else if (tmp == 0)      tmp = 1;

        int i = hydbond_struct.HB_A [tmp];

        ///////////////////////////////////
        /// Loop over all hbonds a particular atom involves
        for (int j = 1; j <= hydbond_struct.HB_ATOM_N [i]; j++) {   //   generates a sequence of bonds for atom "i" 
            int k = hydbond_struct.HB_ATOM [i][j];
            if (hydbond_struct.hb_status [k] == 0) {
                nK += 1;
                hydbond_struct.hb_K [nK] = k;
            }
        }

        ///////////////////////////////////////////////////////
        for (i = 1; i <= nK; i++) {    //  randomly swaps the generated sequence
#ifdef _OPENMP
            r = (double) rand_r (myseed) / RAND_MAX;
#else
            r = (double) rand () / RAND_MAX;
#endif

            double tmp = (int) ceil (r * nK);
            if      (tmp == nK + 1) tmp = nK;
            else if (tmp == 0)      tmp = 1;

            int k = hydbond_struct.hb_K [i];
            hydbond_struct.hb_K [i] = hydbond_struct.hb_K [tmp];
            hydbond_struct.hb_K [tmp] = k;
        }

        //// Choosing which bonds will be turned off
        int k = hydbond_struct.hb_K [1];
        for (i = 2; i <= nK; i++) {
            double j = hydbond_struct.hb_K [i];
            double r2 = exp ((hydbond_struct.hb_energy [j] - hydbond_struct.hb_energy [k]) / T);

#ifdef _OPENMP
            r = (double) rand_r (myseed) / RAND_MAX;
#else
            r = (double) rand () / RAND_MAX;
#endif

            if (r < r2) k = j;
        }

        ///////////////////////////////////////////////////////
        hydbond_struct.hb_status [k] = 1;
        for (int j = 1; j <= hydbond_struct.ATOM_HB_N [k]; j++)
            hydbond_struct.EXCESS [hydbond_struct.ATOM_HB [k][j]] -= 1;
    }
}

////////////////////////////////////////////////////////
void update_force_and_energy (_coord_vel_force_struct &coord_vel_force_struct,
                              _hydbond_struct         &hydbond_struct,
                              _energy_struct          &energy_struct) {
    for (int k = 1; k <= hydbond_struct.hb_N; k++) {
        if (hydbond_struct.hb_status [k] == -1) continue;

        if (hydbond_struct.hb_status [k] == 1) {
            hydbond_struct.hb_status [k] = 0;
            continue;
        }

        if (hydbond_struct.hb_dode [k])
             energy_struct.E3_HB += hydbond_struct.hb_energy [k];
        else energy_struct.E2_HB += hydbond_struct.hb_energy [k];

        double tmp = 0;

        gen_bond_force          (hydbond_struct.hb_k11 [k], hydbond_struct.hb_k21 [k], hydbond_struct.k_r * hydbond_struct.hb_f [k], hydbond_struct.hb_r [k], coord_vel_force_struct, tmp);
        if (not hydbond_struct.ligand_bond [k]) {
            gen_valence_force       (hydbond_struct.hb_k10 [k], hydbond_struct.hb_k11 [k], hydbond_struct.hb_k21 [k], hydbond_struct.k_theta * hydbond_struct.hb_f [k], hydbond_struct.hb_theta1 [k], coord_vel_force_struct, tmp);
            gen_valence_force       (hydbond_struct.hb_k20 [k], hydbond_struct.hb_k21 [k], hydbond_struct.hb_k11 [k], hydbond_struct.k_theta * hydbond_struct.hb_f [k], hydbond_struct.hb_theta2 [k], coord_vel_force_struct, tmp);

            stacking_dihedral_force (hydbond_struct.hb_k10 [k], hydbond_struct.hb_k11 [k], hydbond_struct.hb_k21 [k], hydbond_struct.hb_k20 [k], hydbond_struct.k_phi * hydbond_struct.hb_f [k], hydbond_struct.hb_psi0  [k], coord_vel_force_struct);
            stacking_dihedral_force (hydbond_struct.hb_k21 [k], hydbond_struct.hb_k11 [k], hydbond_struct.hb_k10 [k], hydbond_struct.hb_k12 [k], hydbond_struct.k_phi * hydbond_struct.hb_f [k], hydbond_struct.hb_psi10 [k], coord_vel_force_struct);
            stacking_dihedral_force (hydbond_struct.hb_k11 [k], hydbond_struct.hb_k21 [k], hydbond_struct.hb_k20 [k], hydbond_struct.hb_k22 [k], hydbond_struct.k_phi * hydbond_struct.hb_f [k], hydbond_struct.hb_psi20 [k], coord_vel_force_struct);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
void hydrogen_bonds_interactions (_coord_vel_force_struct &coord_vel_force_struct,
                                  _hydbond_struct         &hydbond_struct,
                                  const double &T,
                                  _energy_struct          &energy_struct,
                                  unsigned int *myseed) {

    compute_hbond (coord_vel_force_struct, hydbond_struct);
    remove_excess_hbond (hydbond_struct, T, myseed);
    update_force_and_energy (coord_vel_force_struct, hydbond_struct, energy_struct);
}

#endif
