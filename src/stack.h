#ifndef STACK_H
#define STACK_H

/////////////////////////////////////////////////////////////
_Dcoordinate norm_base (const _Res_aa &res, const _Dcoordinate &COM) {
    std::string a, b;
    if      (res.name == "A") { a = "C8"; b = "N6"; }
    else if (res.name == "G") { a = "C8"; b = "O6"; }
    else if (res.name == "C") { a = "N4"; b = "O2"; }
    else if (res.name == "U") { a = "O4"; b = "O2"; }

    _Dcoordinate veca, vecb;
    for (size_t i = 0; i < res.atom.size(); i++)
        if (res.atom[i].name == a)
            veca = gen_vector (COM, res.atom[i].coord);
        else if (res.atom[i].name == b)
            vecb = gen_vector (COM, res.atom[i].coord);

    return cross_product (veca, vecb);
}

///////////////////////////////////////////////////////////////////
// Base stacking criteria follows Condon et al, JCTC 2015, 11, 2729
///////////////////////////////////////////////////////////////////
void search_stacking (const _Mol_aa &mol_aa,
                      const _Mol_cg &mol_cg,
//                      const std::vector<int> &ignore_stack_res,
                      _topol_struct &topol_struct,
                      double score_cutoff) {

    for (size_t i = 0; i < mol_cg.res.size(); i++) {
        _Bead bead1 = mol_cg.res[i].bead[2];
        if (bead1.name.empty()) continue;
//        if (bead1.name.empty() || std::find (ignore_stack_res.begin(), ignore_stack_res.end(), mol_cg.res[i].index) != ignore_stack_res.end()) continue;
        _Res_aa res1 = mol_aa.res[i];

        _Dcoordinate norm_base1 = norm_base (res1, bead1.coord);
        for (size_t j = i + 1; j < mol_cg.res.size(); j++) {
            _Bead bead2 = mol_cg.res[j].bead[2];
            if (bead2.name.empty ()) continue;
//            if (bead2.name.empty () || std::find (ignore_stack_res.begin(), ignore_stack_res.end(), mol_cg.res[j].index) != ignore_stack_res.end()) continue;

            bool secondary_stack = 0;
            if (i + 1 == j && mol_cg.res[i].index != topol_struct.ligand_index && mol_cg.res[j].index != topol_struct.ligand_index && \
                not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), res1.index))
                secondary_stack = 1;

            double stack_score = 0;

            ////// distance
            double d = dist_coord (bead1.coord, bead2.coord);
            if      (d <= 3.5) stack_score += 1.;
            else if (d <  5.0) stack_score += 28.583333 * (5.0 - d) / (d*d*d);
            else continue;

            /////// omega
            _Dcoordinate COM1_COM2 = gen_vector (bead1.coord, bead2.coord);
            double omega = compute_angle (norm_base1, COM1_COM2);
            if (omega > 90.) omega = 180. - omega;

            if      (omega <= 25.) stack_score += 1.;
            else if (omega <= 50.) stack_score += 2. - 0.04*omega;
            else continue;

            if (stack_score < score_cutoff) continue; // loosen the criteria a bit (original paper is 1.0)

            ///// xi
            _Res_aa res2 = mol_aa.res[j];
            _Dcoordinate norm_base2 = norm_base (res2, bead2.coord);
            double xi = compute_angle (norm_base1, norm_base2);
            if (xi > 45. && xi < 135.) continue;  // T-shape

            //// Two bases now considered stacked, find out which side of those bases stacks
            int base1 = res1.index;
            if (norm_base1.x*COM1_COM2.x + norm_base1.y*COM1_COM2.y + norm_base1.z*COM1_COM2.z > 0)
                base1 = -base1;
            int base2 = res2.index;
            if (norm_base2.x*COM1_COM2.x + norm_base2.y*COM1_COM2.y + norm_base2.z*COM1_COM2.z < 0)
                base2 = -base2;
            std::pair<int, int> stacked_base = std::make_pair (base1, base2);
//            topol_struct.total_stack.push_back (stacked_base);
            if (secondary_stack) topol_struct.secondary_stack.push_back (stacked_base);
            else                 topol_struct.tertiary_stack.push_back  (stacked_base);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Find which side of base (including the sign defined above) that has only 3rd stacks (and no 2nd stack)
// Later in the simulation, the 3rd stacks have priority over the 2nd stacks
// That is, once the 3rd stacks have formed, do not allow 2nd stacks anymore for this base side
// total_stack stores pairs of stacking bases (eg, -1:::5, 2:::-7)
// stack_status stores pairs of <+/-res_id:::status> where status is 1 if this side only participate in 3rd stackings, and 0 if 2nd stacking is involved
///////////////////////////////////////////////////////////////////////////////
/*void generate_excess_stack (std::vector< std::pair<int, int> >  &total_stack,
                            std::vector< std::pair<int, bool> > &stack_status) {
    for (size_t i = 0; i < total_stack.size(); i++) {
        int id1 = total_stack[i].first;
        if (id1 < 0) id1 = -id1;
        int id2 = total_stack[i].second;
        if (id2 < 0) id2 = -id2;

        /// Secondary stacks
        if (id1 + 1 == id2) {
            // Assign both bases 0 unless they have already existed in the database
            bool found1 = 0, found2 = 0;
            for (size_t j = 0; j < stack_status.size(); j++) {
                if (stack_status[j].first == total_stack[i].first) {
                    stack_status[j].second = 0;
                    found1 = 1;
                }
                if (stack_status[j].first == total_stack[i].second) {
                    stack_status[j].second = 0;
                    found2 = 1;
                }
                if (found1 && found2) break;
            }
            if (not found1) stack_status.push_back (std::make_pair (total_stack[i].first,  0));
            if (not found2) stack_status.push_back (std::make_pair (total_stack[i].second, 0));

        /// Tertiary stacks
        } else {
            // Search thru the database, if not already exist -> assign 1
            bool found1 = 0, found2 = 0;
            for (size_t j = 0; j < stack_status.size(); j++) {
                if (stack_status[j].first == total_stack[i].first)  found1 = 1;
                if (stack_status[j].first == total_stack[i].second) found2 = 1;
                if (found1 && found2) break;
            }
            if (not found1) stack_status.push_back (std::make_pair (total_stack[i].first,  1));
            if (not found2) stack_status.push_back (std::make_pair (total_stack[i].second, 1));
        }
    }
}*/

//////////////////////////////////////////////////////////////////////////////
void initialize_unprocessed_stacks (char * file_unprocessed_stacks,
                                    const _Mol_aa &mol_aa,
                                    const _Mol_cg &mol_cg,
                             //       const std::vector<int> &ignore_stack_res,
                                    _topol_struct &topol_struct,
                                    _stack_struct &stack_struct,
                                    const _coord_vel_force_struct &coord_vel_force_struct,
                                    double T) {    // temperature

    stack_struct.k_r2   = 1.4;         stack_struct.k_r3     = 5.0;
                                       stack_struct.k_theta3 = 1.5;
    stack_struct.k_psi2 = 4.0;         stack_struct.k_psi3   = 0.15;

    search_stacking (mol_aa, mol_cg, topol_struct, stack_struct.score_cutoff);

    ///// Currently not using tertiary stacks detected by the program
    topol_struct.tertiary_stack.clear();

//    stack_struct.ST_ATOM_N = (int *)  calloc (topol_struct.Nres_biomol + 1, sizeof(int));
    stack_struct.ST_ATOM_N.resize (2*topol_struct.Nres_biomol + 1, 0);
//    stack_struct.ST_ATOM   = (int **) calloc (topol_struct.Nres_biomol + 1, sizeof(int *));
//    stack_struct.ST_ATOM.resize   (topol_struct.Nres_biomol + 1);

//    for (int i = 0; i < stack_struct.ST_ATOM.size(); i++)
//        stack_struct.ST_ATOM [i].push_back (0);    // Not using 0 indices

    /////////////// Not using 0 indices
    stack_struct.st_i3.push_back     (0);
    stack_struct.st_j3.push_back     (0);

    stack_struct.st_r.push_back      (0);
    stack_struct.st_theta1.push_back (0);
    stack_struct.st_theta2.push_back (0);
    stack_struct.st_psi.push_back    (0);
    stack_struct.st_psi1.push_back   (0);
    stack_struct.st_psi2.push_back   (0);
    stack_struct.st_E.push_back      (0);
    stack_struct.st_status.push_back (-1);
    stack_struct.ligand_stack.push_back (0);
    ////////////////////////////////

    if (file_unprocessed_stacks != NULL) {
        FILE * f1;
        f1 = fopen (file_unprocessed_stacks, "r");
        fseek (f1, 0L, SEEK_END);
        int finish = ftell (f1);

        stack_struct.s3_N = 0;
        fseek (f1, 0L, SEEK_SET);
        int position;
        do {
            char A1 [10], A2 [10];
            fscanf (f1, "%s %s\n", A1, A2);
            int RES1 = atoi (A1);
            int RES2 = atoi (A2);

            stack_struct.s3_N ++;
            topol_struct.tertiary_stack.push_back (std::make_pair (RES1, RES2));

            // stack_struct.s3_N = topol_struct.tertiary_stack.size();
            /////////// Tertiary stacks
            // for (int index = 0; index < stack_struct.s3_N; index++) {
            //     int RES1 = topol_struct.tertiary_stack[index].first;

            if (RES1 < 0) RES1 = -RES1;
            int res1 = topol_struct.res_pdb2res_simu [RES1];

            // int RES2 = topol_struct.tertiary_stack[index].second;
            if (RES2 < 0) RES2 = -RES2;
            int res2 = topol_struct.res_pdb2res_simu [RES2];

            // int k10 = return_atom_key (topol_struct.phosphate_first, RES1,     topol_struct, "RI");
            // int k11 = return_atom_key (topol_struct.phosphate_first, RES1,     topol_struct, "Base");
            // int k12 = return_atom_key (topol_struct.phosphate_first, RES1 + 1, topol_struct, "PH");
            // int k20 = return_atom_key (topol_struct.phosphate_first, RES2,     topol_struct, "RI");
            // int k21 = return_atom_key (topol_struct.phosphate_first, RES2,     topol_struct, "Base");
            // int k22 = return_atom_key (topol_struct.phosphate_first, RES2 + 1, topol_struct, "PH");

            bool ligand_stack = 0;
            int k10 = topol_struct.atom_key [res1][1];
            int k11 = topol_struct.atom_key [res1][2];
            int k12;
            if (RES1 == topol_struct.ligand_index) {
                k12 = 0;
                ligand_stack = 1;
            } else if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES1) && \
                    res1 < topol_struct.Nres_biomol)
                k12 = topol_struct.atom_key [res1+1][0];
            else {
                printf ("Ignore tertiary stacking between  %d - %d!!! Check the structure !!!\n", RES1, RES2);
                exit (0);
                // if (std::find (topol_struct.stack_ignore.begin(), topol_struct.stack_ignore.end(), RES1) == topol_struct.stack_ignore.end())
                //     topol_struct.stack_ignore.push_back (RES1);
                // continue;
            }

            int k20 = topol_struct.atom_key [res2][1];
            int k21 = topol_struct.atom_key [res2][2];
            int k22;
            if (RES2 == topol_struct.ligand_index) {
                k22 = 0;
                ligand_stack = 1;
            } else if (not std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES2) && \
                    res2 < topol_struct.Nres_biomol)
                k22 = topol_struct.atom_key [res2+1][0];
            else {
                printf ("Ignore tertiary stacking between  %d - %d!!! Check the structure !!!\n", RES1, RES2);
                exit (0);
                // if (std::find (topol_struct.stack_ignore.begin(), topol_struct.stack_ignore.end(), RES2) == topol_struct.stack_ignore.end())
                //     topol_struct.stack_ignore.push_back (RES2);
                // continue;
            }

            stack_struct.ligand_stack.push_back (ligand_stack);

            //////////// just to be able to record "-" and "+" sites separately /////////
            int i = res1 + res1;
            int j = res2 + res2;

            // if (topol_struct.tertiary_stack[index].first  < 0) i--;
            // if (topol_struct.tertiary_stack[index].second < 0) j--;
            if (topol_struct.tertiary_stack.back().first  < 0) i--;
            if (topol_struct.tertiary_stack.back().second < 0) j--;
            stack_struct.st_i3.push_back     (i);
            stack_struct.st_j3.push_back     (j);
            stack_struct.st_r.push_back      (CC_distance    (k11, k21, coord_vel_force_struct));

            if (not ligand_stack) {
                stack_struct.st_theta1.push_back (valence_angle  (k10, k11, k21, coord_vel_force_struct));
                stack_struct.st_theta2.push_back (valence_angle  (k20, k21, k11, coord_vel_force_struct));
                stack_struct.st_psi.push_back    (dihedral_angle (k10, k11, k21, k20, coord_vel_force_struct));
                stack_struct.st_psi1.push_back   (dihedral_angle (k21, k11, k10, k12, coord_vel_force_struct));
                stack_struct.st_psi2.push_back   (dihedral_angle (k11, k21, k20, k22, coord_vel_force_struct));
            } else {
                stack_struct.st_theta1.push_back (0);
                stack_struct.st_theta2.push_back (0);
                stack_struct.st_psi.push_back    (0);
                stack_struct.st_psi1.push_back   (0);
                stack_struct.st_psi2.push_back   (0);
            }

            stack_struct.st_E.push_back      (- stack_struct.st_D);
            stack_struct.st_status.push_back (1);  // Native stack

            stack_struct.ST_ATOM_N [i] ++;
            //        stack_struct.ST_ATOM [i] = (int *) realloc (stack_struct.ST_ATOM [i], (stack_struct.ST_ATOM_N [i] + 1) * sizeof(int));
            //        stack_struct.ST_ATOM [i][stack_struct.ST_ATOM_N [i]] = index+1;

            stack_struct.ST_ATOM_N [j] ++;
            //        stack_struct.ST_ATOM [j] = (int *) realloc (stack_struct.ST_ATOM [j], (stack_struct.ST_ATOM_N [j] + 1) * sizeof(int));
            //        stack_struct.ST_ATOM [j][stack_struct.ST_ATOM_N [j]] = index+1;
            position = ftell (f1);
        } while (position < finish);

        fclose (f1);
    }

    //////////////// Secondary stacks
    stack_struct.st_N = stack_struct.s3_N;  // st_N stores total numer of stacking, starts here and updates later
    stack_struct.st_i2.push_back (0);
    stack_struct.st_j2.push_back (0);

    for (int i = 2; i <= 2*topol_struct.Nres_biomol - 2; i += 2) {
        int res = i / 2;
        int RES = topol_struct.res_simu2res_pdb [res];
        bool next = 0;
        if (topol_struct.atom_key [res][0] == 0 || RES == topol_struct.ligand_index || \
            std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES  )) {
            next = 1;
            if (std::find (topol_struct.stack_ignore.begin(), topol_struct.stack_ignore.end(), RES) == topol_struct.stack_ignore.end())
                topol_struct.stack_ignore.push_back (RES);
        } else {
            int RES2 = topol_struct.res_simu2res_pdb [res + 1];
            if (res + 1 == topol_struct.Nres_biomol || topol_struct.atom_key [res+2][0] == 0 || RES2 == topol_struct.ligand_index || \
                std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES2)) {
                next = 1;
                if (std::find (topol_struct.stack_ignore.begin(), topol_struct.stack_ignore.end(), RES2) == topol_struct.stack_ignore.end())
                    topol_struct.stack_ignore.push_back (RES2);
            }
        }

        if (next) continue;

        stack_struct.st_N ++;
        stack_struct.st_i2.push_back (i);
        stack_struct.st_j2.push_back (i+1);

        ///////// Checking for native stack ///////////
        int nextRES = -RES-1;
        bool native = 0;
        for (size_t j = 0; j < topol_struct.secondary_stack.size(); j++)
            if ((topol_struct.secondary_stack[j].first == RES     && topol_struct.secondary_stack[j].second == nextRES) || \
                (topol_struct.secondary_stack[j].first == nextRES && topol_struct.secondary_stack[j].second == RES)) {
                native = 1;
                break;
            }
        stack_struct.st_status.push_back (native);

        ////////// Check those HARD-CODE
        if (mol_cg.res[res - 1].name == "A") {
            if (mol_cg.res[res].name == "A") {
                stack_struct.st_r.push_back (4.1806530);
                stack_struct.st_E.push_back (-5.194 + stack_struct.ss_D / 0.7093838769 - 0.319*(T - 0.594));

            } else if (mol_cg.res[res].name == "C") {
                stack_struct.st_r.push_back (3.8260185);
                stack_struct.st_E.push_back (-5.146 + stack_struct.ss_D / 0.7185955600 - 0.319*(T - 0.594));

            } else if (mol_cg.res[res].name == "G") {
                stack_struct.st_r.push_back (4.4255305);
                stack_struct.st_E.push_back (-5.977 + stack_struct.ss_D / 0.6968019829 + 5.301*(T - 0.678));

            } else if (mol_cg.res[res].name == "U") {
                stack_struct.st_r.push_back (3.8260185);
                stack_struct.st_E.push_back (-5.146 + stack_struct.ss_D / 0.7185955600 - 0.319*(T - 0.594));
            }

        } else if (mol_cg.res[res - 1].name == "C") {
            if (mol_cg.res[res].name == "A") {
                stack_struct.st_r.push_back (4.7010580);
                stack_struct.st_E.push_back (-5.163 + stack_struct.ss_D / 0.6847830171 - 0.319*(T - 0.594));

            } else if (mol_cg.res[res].name == "C") {
                stack_struct.st_r.push_back (4.2500910);
                stack_struct.st_E.push_back (-4.873 + stack_struct.ss_D / 0.6991615586 - 1.567*(T - 0.568));

            } else if (mol_cg.res[res].name == "G") {
                stack_struct.st_r.push_back (4.9790760);
                stack_struct.st_E.push_back (-5.482 + stack_struct.ss_D / 0.6816268897 + 0.774*(T - 0.627));

            } else if (mol_cg.res[res].name == "U") {
                stack_struct.st_r.push_back (4.2273615);
                stack_struct.st_E.push_back (-4.873 + stack_struct.ss_D / 0.6832570771 - 1.567*(T - 0.568));
            }

        } else if (mol_cg.res[res - 1].name == "G") {
            if (mol_cg.res[res].name == "A") {
                stack_struct.st_r.push_back (4.0128560);
                stack_struct.st_E.push_back (-5.948 + stack_struct.ss_D / 0.6903176657 + 5.301*(T - 0.678));

            } else if (mol_cg.res[res].name == "C") {
                stack_struct.st_r.push_back (3.6784360);
                stack_struct.st_E.push_back (-5.927 + stack_struct.ss_D / 0.7042060343 + 4.370*(T - 0.682));

            } else if (mol_cg.res[res].name == "G") {
                stack_struct.st_r.push_back (4.2427250);
                stack_struct.st_E.push_back (-6.416 + stack_struct.ss_D / 0.6971421514 + 7.346*(T - 0.728));

            } else if (mol_cg.res[res].name == "U") {
                stack_struct.st_r.push_back (3.6616930);
                stack_struct.st_E.push_back (-5.836 + stack_struct.ss_D / 0.6984426543 + 2.924*(T - 0.672));
            }

        } else if (mol_cg.res[res - 1].name == "U") {
            if (mol_cg.res[res].name == "A") {
                stack_struct.st_r.push_back (4.7010580);
                stack_struct.st_E.push_back (-5.163 + stack_struct.ss_D / 0.6847830171 - 0.319*(T - 0.594));

            } else if (mol_cg.res[res].name == "C") {
                stack_struct.st_r.push_back (4.2679180);
                stack_struct.st_E.push_back (-4.880 + stack_struct.ss_D / 0.6758595771 - 1.567*(T - 0.568));

            } else if (mol_cg.res[res].name == "G") {
                stack_struct.st_r.push_back (4.9977560);
                stack_struct.st_E.push_back (-5.886 + stack_struct.ss_D / 0.7025528229 + 2.924*(T - 0.672));

            } else if (mol_cg.res[res].name == "U") {
                stack_struct.st_r.push_back (4.2453650);
                stack_struct.st_E.push_back (-4.267 + stack_struct.ss_D / 0.6686014771 - 3.563*(T - 0.500));
            }
        }

        stack_struct.ST_ATOM_N [i] ++;
//        stack_struct.ST_ATOM [i] = (int *) realloc (stack_struct.ST_ATOM [i], (stack_struct.ST_ATOM_N [i] + 1) * sizeof(int));
//        stack_struct.ST_ATOM [i][stack_struct.ST_ATOM_N [i]] = stack_struct.st_N;
//        stack_struct.ST_ATOM [i].push_back (stack_struct.st_N);

        stack_struct.ST_ATOM_N [i + 1] ++;
//        stack_struct.ST_ATOM [i + 1] = (int *) realloc (stack_struct.ST_ATOM [i + 1], (stack_struct.ST_ATOM_N [i + 1] + 1) * sizeof(int));
//        stack_struct.ST_ATOM [i + 1][stack_struct.ST_ATOM_N [i + 1]] = stack_struct.st_N;
//        stack_struct.ST_ATOM [i + 1].push_back (stack_struct.st_N);
    }

    if (topol_struct.stack_ignore.size() > 0) {
        printf ("Ignoring stacking from those nucleotides: ");
        for (size_t i = 0; i < topol_struct.stack_ignore.size(); i++)
            printf ("%d  ", topol_struct.stack_ignore [i]);
        printf ("\n");
    }

    stack_struct.st_psi0.resize   (stack_struct.s3_N + 1, 0);
    stack_struct.st_psi10.resize  (stack_struct.st_N + 1, 0);
    stack_struct.st_psi20.resize  (stack_struct.st_N + 1, 0);
    stack_struct.st_energy.resize (stack_struct.st_N + 1, 0);
    stack_struct.st_status.resize (stack_struct.st_N + 1, 1);

    stack_struct.e3_stack.resize  (stack_struct.s3_N + 1, 0);
//    stack_struct.ST_EXCESS.resize (topol_struct.Nres_biomol + 1, 0);
//    stack_struct.ST_EXCESS.resize (2*topol_struct.Nres_biomol + 1, 0);
    stack_struct.Nstack_allowed.resize (2*topol_struct.Nres_biomol + 1, 0);
    stack_struct.Nstack_formed .resize (2*topol_struct.Nres_biomol + 1, 0);

    // Set the maximum stacking interactions one side of base can participate
    std::vector<std::pair<int, int> > total_stack = topol_struct.secondary_stack;
    total_stack.insert (total_stack.end(), topol_struct.tertiary_stack.begin(), topol_struct.tertiary_stack.end());
//    for (int i = 0; i < total_stack.size(); i++)
//        printf ("RES1 %d   - RES2 %d\n", total_stack [i].first, total_stack [i].second);

    for (size_t i = 0; i < total_stack.size(); i++) {
        int RES1 = total_stack [i].first;
        if (RES1 < 0) RES1 = -RES1;
        int res1 = topol_struct.res_pdb2res_simu [RES1];
        int stack_index = res1 + res1;
        if (total_stack [i].first < 0) stack_index--;
        stack_struct.Nstack_allowed [stack_index] ++;

        int RES2 = total_stack [i].second;
        if (RES2 < 0) RES2 = -RES2;
        int res2 = topol_struct.res_pdb2res_simu [RES2];
        stack_index = res2 + res2;
        if (total_stack [i].second < 0) stack_index--;
        stack_struct.Nstack_allowed [stack_index] ++;
    }

    // Allow at least one stacking per base-side
    for (int i = 1; i <= 2*topol_struct.Nres_biomol; i++)
        if (stack_struct.Nstack_allowed [i] == 0)
            stack_struct.Nstack_allowed [i] = 1;

    printf ("Double-stacking base: ");
    for (int i = 1; i <= 2*topol_struct.Nres_biomol; i++)
        if (stack_struct.Nstack_allowed [i] == 2) {
            int res = (i / 2) + (i % 2 == 1);
            int RES = topol_struct.res_simu2res_pdb [res];
            if (i % 2 == 1) RES = -RES;
            printf ("%d  ", RES);
        }
    printf ("\n");
}

/////////////////////////////////////////////////////////////////
void stacking_interactions (const _topol_struct     &topol_struct,
                            _coord_vel_force_struct &coord_vel_force_struct,
                            _stack_struct           &stack_struct,
                            _energy_struct          &energy_struct) {
    //// Start at zero for all sides
    std::fill (stack_struct.Nstack_formed.begin(), stack_struct.Nstack_formed.end(), 0);

    ////// Tertiary stacks
    for (int k = 1; k <= stack_struct.s3_N; k++) {
        int i = stack_struct.st_i3 [k];
        int j = stack_struct.st_j3 [k];

        int res1 = (i / 2) + (i % 2 == 1);
        int res2 = (j / 2) + (j % 2 == 1);

//        int k10 = topol_struct.atom_key [i][0];
//        int k11 = topol_struct.atom_key [i][1];
//        int k12 = topol_struct.atom_key [i + 1][0];
//        int k20 = topol_struct.atom_key [j][0];
//        int k21 = topol_struct.atom_key [j][1];
//        int k22 = topol_struct.atom_key [j + 1][0];

        int k10 = topol_struct.atom_key [res1]  [1];
        int k11 = topol_struct.atom_key [res1]  [2];
        int k12, k22;

        int k20 = topol_struct.atom_key [res2]  [1];
        int k21 = topol_struct.atom_key [res2]  [2];

        if (not stack_struct.ligand_stack [k]) {
            k12 = topol_struct.atom_key [res1+1][0];
            k22 = topol_struct.atom_key [res2+1][0];
        } else {
            k12 = 0;
            k22 = 0;
        }

        double r = CC_distance (k11, k21, coord_vel_force_struct);
        stack_struct.e3_stack [k] = r - stack_struct.st_r [k];

        // Once tertiary stacks form, do not consider secondary stacks
        if (fabs(stack_struct.e3_stack [k]) < 10.0) {
            stack_struct.Nstack_formed [i]++;
            stack_struct.Nstack_formed [j]++;
        }

        double r2 = 1. + stack_struct.k_r3 * stack_struct.e3_stack [k] * stack_struct.e3_stack [k];

        if (not stack_struct.ligand_stack [k]) {
            double theta1 = valence_angle (k10, k11, k21, coord_vel_force_struct);
            double theta2 = valence_angle (k20, k21, k11, coord_vel_force_struct);

            double psi  = dihedral_angle (k10, k11, k21, k20, coord_vel_force_struct);
            double psi1 = dihedral_angle (k21, k11, k10, k12, coord_vel_force_struct);
            double psi2 = dihedral_angle (k11, k21, k20, k22, coord_vel_force_struct);

            if (fabs (psi - stack_struct.st_psi [k]) <= PI)
                stack_struct.st_psi0 [k] = stack_struct.st_psi [k];
            else if (psi - stack_struct.st_psi [k] > PI)
                stack_struct.st_psi0 [k] = stack_struct.st_psi [k] + 2.*PI;
            else stack_struct.st_psi0 [k] = stack_struct.st_psi [k] - 2.*PI;

            if (fabs (psi1 - stack_struct.st_psi1 [k]) <= PI)
                stack_struct.st_psi10 [k] = stack_struct.st_psi1 [k];
            else if (psi1 - stack_struct.st_psi1 [k] > PI)
                stack_struct.st_psi10 [k] = stack_struct.st_psi1 [k] + 2.*PI;
            else stack_struct.st_psi10 [k] = stack_struct.st_psi1 [k] - 2.*PI;

            if (fabs (psi2 - stack_struct.st_psi2 [k]) <= PI)
                stack_struct.st_psi20 [k] = stack_struct.st_psi2 [k];
            else if (psi2 - stack_struct.st_psi2 [k] > PI)
                stack_struct.st_psi20 [k] = stack_struct.st_psi2 [k] + 2.*PI;
            else stack_struct.st_psi20 [k] = stack_struct.st_psi2 [k] - 2.*PI;

            ///// HARD-CODE
            r2 += stack_struct.k_theta3 * (theta1 - stack_struct.st_theta1 [k]) * (theta1 - stack_struct.st_theta1 [k]);
            r2 += stack_struct.k_theta3 * (theta2 - stack_struct.st_theta2 [k]) * (theta2 - stack_struct.st_theta2 [k]);
            r2 += stack_struct.k_psi3   * (psi  - stack_struct.st_psi0  [k]) * (psi  - stack_struct.st_psi0  [k]);
            r2 += stack_struct.k_psi3   * (psi1 - stack_struct.st_psi10 [k]) * (psi1 - stack_struct.st_psi10 [k]);
            r2 += stack_struct.k_psi3   * (psi2 - stack_struct.st_psi20 [k]) * (psi2 - stack_struct.st_psi20 [k]);
        }
        stack_struct.st_energy [k] = stack_struct.st_E [k] / r2;
        energy_struct.E3_STACK += stack_struct.st_energy [k];

        double st_f = - stack_struct.st_energy [k] / r2;

        double tmp = 0;
        gen_bond_force          (k11, k21,           stack_struct.k_r3 * st_f, stack_struct.st_r [k], coord_vel_force_struct, tmp);
        if (not stack_struct.ligand_stack [k]) {
            gen_valence_force       (k10, k11, k21,      stack_struct.k_theta3 * st_f, stack_struct.st_theta1 [k], coord_vel_force_struct, tmp);
            gen_valence_force       (k20, k21, k11,      stack_struct.k_theta3 * st_f, stack_struct.st_theta2 [k], coord_vel_force_struct, tmp);
            stacking_dihedral_force (k10, k11, k21, k20, stack_struct.k_psi3   * st_f, stack_struct.st_psi0 [k],  coord_vel_force_struct);
            stacking_dihedral_force (k21, k11, k10, k12, stack_struct.k_psi3   * st_f, stack_struct.st_psi10 [k], coord_vel_force_struct);
            stacking_dihedral_force (k11, k21, k20, k22, stack_struct.k_psi3   * st_f, stack_struct.st_psi20 [k], coord_vel_force_struct);
        }
    }

    ////////// Secondary stacks
    // Not having ligand stacks
    for (int k = stack_struct.s3_N + 1; k <= stack_struct.st_N; k++) {
        int i = stack_struct.st_i2 [k - stack_struct.s3_N];
        int res = (i / 2) + (i % 2 == 1);
        if (stack_struct.Nstack_formed [i]   == stack_struct.Nstack_allowed [i] || \
            stack_struct.Nstack_formed [i+1] == stack_struct.Nstack_allowed [i+1]) continue;

//        int k_10 = topol_struct.atom_key [i - 1][0];
//        int k00  = topol_struct.atom_key [i][0];
//        int k01  = topol_struct.atom_key [i][1];
//        int k10  = topol_struct.atom_key [i + 1][0];
//        int k20  = topol_struct.atom_key [i + 2][0];
//        int k21  = topol_struct.atom_key [i + 2][1];
//        int k30  = topol_struct.atom_key [i + 3][0];

        int k_10 = topol_struct.atom_key [res]  [0];
        int k00  = topol_struct.atom_key [res]  [1];
        int k01  = topol_struct.atom_key [res]  [2];
        int k10  = topol_struct.atom_key [res+1][0];
        int k20  = topol_struct.atom_key [res+1][1];
        int k21  = topol_struct.atom_key [res+1][2];
        int k30  = topol_struct.atom_key [res+2][0];

        double r = CC_distance (k01, k21, coord_vel_force_struct);
        double psi1 = dihedral_angle (k_10, k00, k10, k20, coord_vel_force_struct);
        double psi2 = dihedral_angle (k30, k20, k10, k00, coord_vel_force_struct);

        //////// Check those HARD-CODE
        if (fabs (psi1 + 2.58684) <= PI)
            stack_struct.st_psi10 [k] = - 2.58684;
        else if (psi1 + 2.58684 > PI)
            stack_struct.st_psi10 [k] = - 2.58684 + 2*PI;
        else stack_struct.st_psi10 [k] = - 2.58684 - 2*PI;

        if (fabs (psi2 - 3.07135) <= PI)
            stack_struct.st_psi20 [k] = 3.07135;
        else if (psi2 - 3.07135 > PI)
            stack_struct.st_psi20 [k] = 3.07135 + 2*PI;
        else stack_struct.st_psi20 [k] = 3.07135 - 2*PI;

        double r2 = 1. + stack_struct.k_r2 * (r - stack_struct.st_r [k]) * (r - stack_struct.st_r [k]);
        r2 += stack_struct.k_psi2 * (psi1 - stack_struct.st_psi10 [k]) * (psi1 - stack_struct.st_psi10 [k]);
        r2 += stack_struct.k_psi2 * (psi2 - stack_struct.st_psi20 [k]) * (psi2 - stack_struct.st_psi20 [k]);

        stack_struct.st_energy [k] = stack_struct.st_E [k] / r2;

        double st_f = - stack_struct.st_energy [k] / r2;

        double tmp = 0;
        gen_bond_force          (k01, k21,            stack_struct.k_r2 * st_f, stack_struct.st_r [k], coord_vel_force_struct, tmp);
        stacking_dihedral_force (k_10, k00, k10, k20, stack_struct.k_psi2 * st_f, stack_struct.st_psi10 [k], coord_vel_force_struct);
        stacking_dihedral_force (k30, k20, k10, k00,  stack_struct.k_psi2 * st_f, stack_struct.st_psi20 [k], coord_vel_force_struct);

        if (stack_struct.st_status [k]) energy_struct.E3_STACK += stack_struct.st_energy [k];
        else                            energy_struct.E2_STACK += stack_struct.st_energy [k];
    }
}

#endif
