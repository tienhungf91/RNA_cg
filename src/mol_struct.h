#ifndef MOL_STRUCT_H
#define MOL_STRUCT_H

/////////////////////////////////////////////////////////////////////
/*_Mol_cg SOP_Coarse_graining (const _Mol_aa &mol_aa) {
    _Mol_cg result;
    for (size_t i = 0; i < mol_aa.res.size(); i++) {
        _Res_aa res = mol_aa.res[i];
        _Bead bead;
        bead.connectedTo = res.connectedTo;
        bead.name = res.name;
        bead.index = res.index;
        for (size_t j = 0; j < res.atom.size(); j++) {
            bead.mass += res.atom[j].mass;
            if (res.atom[j].name == "CA")
                bead.coord = res.atom[j].coord;
        }
        bead.name = "C";
        result.bead.push_back (bead);
    }
    result.type = mol_aa.type;
    return result;
}*/

///////////////////////////////////////////////////////////////////
_Mol_cg TIS_Coarse_graining (const _Mol_aa &mol_aa,
                             const std::vector<int> &ter_card,
                             int ligand_index) {
                         //    _topol_struct &topol_struct) {
    _Mol_cg result;
    result.netcharge = 0;
    size_t bead_index = 0;

    for (size_t i = 0; i < mol_aa.res.size(); i++) {
        _Res_aa res = mol_aa.res[i];
        _Res_cg res_cg;
        res_cg.name = res.name;
        res_cg.index = res.index;
        res_cg.chain = res.chain;
        res_cg.connectedTo = res.connectedTo;

        _Bead phosphate, ribose, base;
        res_cg.bead.resize (3, phosphate); // initialize null beads

        ribose.name = "RI";
        ribose.coord.x = 0;   ribose.coord.y = 0;   ribose.coord.z = 0;

        base.name = res.name;
        base.coord.x = 0;     base.coord.y = 0;     base.coord.z = 0;

        phosphate.name = "PH";
        phosphate.coord.x = 0;   phosphate.coord.y = 0;   phosphate.coord.z = 0;

        size_t count_R = 0, count_B = 0, count_P = 0;
        for (size_t j = 0; j < res.atom.size(); j++) {
            std::string atm_name = res.atom[j].name;

            // Center of geometry
            if (atm_name.substr(0, 1) != "H") {
                if (atm_name == "C1'" || atm_name == "C2'" || atm_name == "C3'" || \
                    atm_name == "C4'" || atm_name == "C5'" || atm_name == "O2'" || \
                    atm_name == "O3'" || atm_name == "O4'" || atm_name == "O5'") {
                    count_R++;
                    ribose.coord.x += res.atom[j].coord.x;
                    ribose.coord.y += res.atom[j].coord.y;
                    ribose.coord.z += res.atom[j].coord.z;
                } else if (atm_name == "P" || atm_name == "OP1" || atm_name == "OP2" || \
                                              atm_name == "O1P" || atm_name == "O2P") {
                    count_P++;
                    phosphate.coord.x += res.atom[j].coord.x;
                    phosphate.coord.y += res.atom[j].coord.y;
                    phosphate.coord.z += res.atom[j].coord.z;
                } else {
                    count_B++;
                    base.coord.x += res.atom[j].coord.x;
                    base.coord.y += res.atom[j].coord.y;
                    base.coord.z += res.atom[j].coord.z;
                }
            }
        }

        if (count_P == 3) {
            phosphate.coord.x /= count_P;
            phosphate.coord.y /= count_P;
            phosphate.coord.z /= count_P;

            bead_index++;
            phosphate.index = bead_index;
            res_cg.charge = -1;
            result.netcharge--;
        } else if (count_P > 0) {
            printf ("Nucleotide %d not having a complete phosphate group!!\n", res.index);
            exit (0);
        }

        if (count_R == 9 && mol_aa.type == "RNA" || count_R == 8 && mol_aa.type == "DNA") {
            ribose.coord.x /= count_R;
            ribose.coord.y /= count_R;
            ribose.coord.z /= count_R;

            bead_index++;
            ribose.index = bead_index;
        } else if (count_R > 0) {
            printf ("Nucleotide %d not having a complete ribose group!!\n", res.index);
            exit (0);
        }

        if ((count_B == 11 && (res.name == "G" || res.name == "DG")) || \
            (count_B ==  8 && (res.name == "C" || res.name == "DC"   || res.name == "U")) || \
            (count_B == 10 && (res.name == "A" || res.name == "DA")) || \
            (count_B ==  9 && res.name == "DT")) {
            base.coord.x /= count_B;
            base.coord.y /= count_B;
            base.coord.z /= count_B;

            bead_index++;
            base.index = bead_index;
        } else if (count_B > 0) {
            printf ("Nucleotide %d not having a complete base group!!\n", res.index);
            exit (0);
        }

        //////////// Connect
        if (count_P == 3 && (count_R == 9 || count_R == 8)) {
            ribose.connectedTo.push_back    (phosphate.index);
            phosphate.connectedTo.push_back (ribose.index);
        }

        if ((count_R == 9 || count_R == 8) && count_B > 0) {
            base.connectedTo.push_back   (ribose.index);
            ribose.connectedTo.push_back (base.index);
        } else if ((res.index != ligand_index) && (count_R == 0 xor count_B == 0)) {
            printf ("Nucleotide %d not having either ribose or base group!!\n", res.index);
            exit (0);
        }

        if (count_P == 3 && result.res.size() > 0 && \
            not std::binary_search (ter_card.begin(), ter_card.end(), result.res.back().index)) {

            result.res.back().bead[1].connectedTo.push_back (phosphate.index);
            phosphate.connectedTo.push_back (result.res.back().bead[1].index);
        }

        if (count_P == 3) res_cg.bead[0] = phosphate;
        if (count_R == 9 || count_R == 8) res_cg.bead[1] = ribose;
        if (count_B > 0)  res_cg.bead[2] = base;
        result.res.push_back (res_cg);
    }
    result.type = mol_aa.type;
    return result;
}

///////////////////////////////////////////////////////////////
_Mol_cg TIS_Coarse_graining_ion (const _Mol_aa &mol_aa,
                                 const std::vector<int> &ter_card,
                                 int ligand_index) {
                             //    _topol_struct &topol_struct) {
    _Mol_cg result;
    result.netcharge = 0;
    size_t bead_index = 0;

    for (size_t i = 0; i < mol_aa.res.size(); i++) {
        _Res_aa res = mol_aa.res[i];
        _Res_cg res_cg;
        res_cg.name = res.name;
        res_cg.index = res.index;
        res_cg.chain = res.chain;
        res_cg.connectedTo = res.connectedTo;

        _Bead phosphate, ribose, base, mg;
        if (res.name != "MG")
            res_cg.bead.resize (3, phosphate); // initialize null beads
        else
            res_cg.bead.resize (1, phosphate);

        ribose.name = "RI";
        ribose.coord.x = 0;   ribose.coord.y = 0;   ribose.coord.z = 0;

        base.name = res.name;
        base.coord.x = 0;     base.coord.y = 0;     base.coord.z = 0;

        phosphate.name = "PH";
        phosphate.coord.x = 0;   phosphate.coord.y = 0;   phosphate.coord.z = 0;

        mg.name = "MG";
        mg.coord.x = 0;    mg.coord.y = 0;     mg.coord.z = 0;

        size_t count_R = 0, count_B = 0, count_P = 0;
        for (size_t j = 0; j < res.atom.size(); j++) {
            std::string atm_name = res.atom[j].name;
            //if (res.name == "MG")
            //    printf ("%s\n", atm_name.c_str());

            // Center of geometry
            if (atm_name.substr(0, 1) != "H") {
                if (atm_name == "C1'" || atm_name == "C2'" || atm_name == "C3'" || \
                    atm_name == "C4'" || atm_name == "C5'" || atm_name == "O2'" || \
                    atm_name == "O3'" || atm_name == "O4'" || atm_name == "O5'") {
                    count_R++;
                    ribose.coord.x += res.atom[j].coord.x;
                    ribose.coord.y += res.atom[j].coord.y;
                    ribose.coord.z += res.atom[j].coord.z;
                } else if (atm_name == "P" || atm_name == "OP1" || atm_name == "OP2" || \
                                              atm_name == "O1P" || atm_name == "O2P") {
                    count_P++;
                    phosphate.coord.x += res.atom[j].coord.x;
                    phosphate.coord.y += res.atom[j].coord.y;
                    phosphate.coord.z += res.atom[j].coord.z;
                } else if (atm_name == "MG") {
                    mg.coord.x += res.atom[j].coord.x;
                    mg.coord.y += res.atom[j].coord.y;
                    mg.coord.z += res.atom[j].coord.z;
                } else {
                    count_B++;
                    base.coord.x += res.atom[j].coord.x;
                    base.coord.y += res.atom[j].coord.y;
                    base.coord.z += res.atom[j].coord.z;
                }
            }
        }

        if (res.name != "MG") {
            if (count_P == 3) {
                phosphate.coord.x /= count_P;
                phosphate.coord.y /= count_P;
                phosphate.coord.z /= count_P;

                bead_index++;
                phosphate.index = bead_index;
                res_cg.charge = -1;
                result.netcharge--;
            } else if (count_P > 0) {
                printf ("Nucleotide %d not having a complete phosphate group!!\n", res.index);
                exit (0);
            }

            if (count_R == 9 && mol_aa.type == "RNA" || count_R == 8 && mol_aa.type == "DNA") {
                ribose.coord.x /= count_R;
                ribose.coord.y /= count_R;
                ribose.coord.z /= count_R;

                bead_index++;
                ribose.index = bead_index;
            } else if (count_R > 0) {
                printf ("Nucleotide %d not having a complete ribose group!!\n", res.index);
                exit (0);
            }

            if ((count_B == 11 && (res.name == "G" || res.name == "DG")) || \
                (count_B ==  8 && (res.name == "C" || res.name == "DC"   || res.name == "U")) || \
                (count_B == 10 && (res.name == "A" || res.name == "DA")) || \
                (count_B ==  9 && res.name == "DT")) {
                base.coord.x /= count_B;
                base.coord.y /= count_B;
                base.coord.z /= count_B;

                bead_index++;
                base.index = bead_index;
            } else if (count_B > 0) {
                printf ("Nucleotide %d not having a complete base group!!\n", res.index);
                exit (0);
            }

            //////////// Connect
            if (count_P == 3 && (count_R == 9 || count_R == 8)) {
                ribose.connectedTo.push_back    (phosphate.index);
                phosphate.connectedTo.push_back (ribose.index);
            }

            if ((count_R == 9 || count_R == 8) && count_B > 0) {
                base.connectedTo.push_back   (ribose.index);
                ribose.connectedTo.push_back (base.index);
            } else if ((res.index != ligand_index) && (count_R == 0 xor count_B == 0)) {
                printf ("Nucleotide %d not having either ribose or base group!!\n", res.index);
                exit (0);
            }

            if (count_P == 3 && result.res.size() > 0 && \
                not std::binary_search (ter_card.begin(), ter_card.end(), result.res.back().index)) {

                result.res.back().bead[1].connectedTo.push_back (phosphate.index);
                phosphate.connectedTo.push_back (result.res.back().bead[1].index);
            }

            if (count_P == 3) res_cg.bead[0] = phosphate;
            if (count_R == 9 || count_R == 8) res_cg.bead[1] = ribose;
            if (count_B > 0)  res_cg.bead[2] = base;

        } else
            res_cg.bead[0] = mg;
        
        result.res.push_back (res_cg);
    }
    result.type = mol_aa.type;
    return result;
}

/////////////////////////////////////////////////////////////
void Parse_mol_struct_to_simu (const _Mol_cg           &mol_cg,
                               _topol_struct           &topol_struct,
                               const std::vector<int>  &res_zeroQ,
                               _coord_vel_force_struct &coord_vel_force_struct) {

//    topol_struct.amino_key.push_back (' ');
    topol_struct.part_key.push_back  (-1);
    topol_struct.maxi_key.push_back  (0);
    topol_struct.Nsite.push_back     (-1);
    topol_struct.res_simu2res_pdb.push_back (-1);
//    topol_struct.atom_key2res_pdb.push_back (std::make_pair (0, 0));

    _Dcoordinate COM;    COM.x = 0;   COM.y = 0;    COM.z = 0;
    topol_struct.Natm_biomol = 0;

    coord_vel_force_struct.coordx.push_back (0);
    coord_vel_force_struct.coordy.push_back (0);
    coord_vel_force_struct.coordz.push_back (0);

    for (size_t i = 0; i < mol_cg.res.size(); i++) {
        _Res_cg res = mol_cg.res[i];

        topol_struct.Nres_biomol++;
        if (topol_struct.res_pdb2res_simu.size() < res.index + 1)
            topol_struct.res_pdb2res_simu.resize (res.index + 1, 0);
        topol_struct.res_pdb2res_simu [res.index] = topol_struct.Nres_biomol;
        topol_struct.res_simu2res_pdb.push_back (res.index);

        topol_struct.atom_key = (int **) realloc (topol_struct.atom_key, (topol_struct.Nres_biomol + 1) * sizeof(int*));
        topol_struct.atom_key [topol_struct.Nres_biomol] = (int *) calloc (3, sizeof(int));

        for (size_t j = 0; j < 3; j++) {
            _Bead bead = res.bead[j];
            if (bead.name.empty ()) continue;

            topol_struct.Natm_biomol++;
            if (res.index == topol_struct.ligand_index) topol_struct.ligand_bead = topol_struct.Natm_biomol;

            //coord_vel_force_struct.coordx = (double *) realloc (coord_vel_force_struct.coordx, (topol_struct.Natm_biomol + 1) * sizeof(double));
            //coord_vel_force_struct.coordy = (double *) realloc (coord_vel_force_struct.coordy, (topol_struct.Natm_biomol + 1) * sizeof(double));
            //coord_vel_force_struct.coordz = (double *) realloc (coord_vel_force_struct.coordz, (topol_struct.Natm_biomol + 1) * sizeof(double));

            //coord_vel_force_struct.coordx [topol_struct.Natm_biomol] = bead.coord.x;
            //coord_vel_force_struct.coordy [topol_struct.Natm_biomol] = bead.coord.y;
            //coord_vel_force_struct.coordz [topol_struct.Natm_biomol] = bead.coord.z;

            coord_vel_force_struct.coordx.push_back (bead.coord.x);
            coord_vel_force_struct.coordy.push_back (bead.coord.y);
            coord_vel_force_struct.coordz.push_back (bead.coord.z);
            COM.x += bead.coord.x;
            COM.y += bead.coord.y;
            COM.z += bead.coord.z;

            if (bead.name == "PH") {    // phosphate
                topol_struct.atom_key [topol_struct.Nres_biomol][0] = topol_struct.Natm_biomol;
                //topol_struct.atom_key2res_pdb.push_back (std::make_pair (topol_struct.Nres_biomol, 0));
                topol_struct.part_key.push_back (0);
                topol_struct.maxi_key.push_back (1);
            } else if (bead.name == "RI") {   // ribose
                topol_struct.atom_key [topol_struct.Nres_biomol][1] = topol_struct.Natm_biomol;
                //topol_struct.atom_key2res_pdb.push_back (std::make_pair (topol_struct.Nres_biomol, 1));
                topol_struct.part_key.push_back (1);
                topol_struct.maxi_key.push_back (2);
            } else {           // base
                topol_struct.atom_key [topol_struct.Nres_biomol][2] = topol_struct.Natm_biomol;
                //topol_struct.atom_key2res_pdb.push_back (std::make_pair (topol_struct.Nres_biomol, 2));
                topol_struct.part_key.push_back (2);
                if      (bead.name == "A" || bead.name == "DA") topol_struct.maxi_key.push_back (3);
                else if (bead.name == "G" || bead.name == "DG") topol_struct.maxi_key.push_back (4);
                else if (bead.name == "C" || bead.name == "DC") topol_struct.maxi_key.push_back (5);
                //// TODO: find parameter for T
                else if (bead.name == "U" || bead.name == "DT") topol_struct.maxi_key.push_back (6);
            }
        }
    }

    for (int i = 0; i < res_zeroQ.size(); i++) {
        int internal_index = topol_struct.res_pdb2res_simu [res_zeroQ[i]];
        topol_struct.bead_zeroQ.push_back (topol_struct.atom_key [internal_index][0]);
    }

    COM.x /= topol_struct.Natm_biomol;
    COM.y /= topol_struct.Natm_biomol;
    COM.z /= topol_struct.Natm_biomol;

    // Find out which bead is close to COM
    double min = 9e9;
    for (size_t i = 0; i < mol_cg.res.size(); i++)
        for (size_t j = 0; j < 3; j++) {
            if (mol_cg.res[i].bead[j].name.empty()) continue;
            _Dcoordinate dist;
            dist.x = COM.x - mol_cg.res[i].bead[j].coord.x;
            dist.y = COM.y - mol_cg.res[i].bead[j].coord.y;
            dist.z = COM.z - mol_cg.res[i].bead[j].coord.z;
            double d = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
            if (d < min) {
                topol_struct.bead_center = topol_struct.atom_key [i+1][j];
                min = d;
            }
        }

/*    topol_struct.INDX.resize (topol_struct.Natm_biomol + 1, -1);
    topol_struct.JNDX.resize (topol_struct.Natm_biomol + 1, -1);

    for (int i = 1; i <= topol_struct.Nres_biomol; i++)
        for (int j = 0; j <= topol_struct.Nsite [i]; j++) {
            int k = topol_struct.atom_key [i][j];
            topol_struct.INDX [k] = i;
            topol_struct.JNDX [k] = j;
        }*/
}

////////////////////////////////////////////////////////////////////////////
// Two beads are considered connected if linked by either 0 or 1 other bead
//     (bond and angle restraint ...)
////////////////////////////////////////////////////////////////////////////
void gen_con_matrix (_topol_struct &topol_struct,
                     const _Mol_cg &mol_cg) {

    //printf ("Number of atoms  %d\n", topol_struct.Natm);
    topol_struct.con_matrix = (int**) calloc (topol_struct.Natm + 1, sizeof(int*));
    for (int i = 0; i <= topol_struct.Natm; i++)
        topol_struct.con_matrix[i] = (int*) calloc (topol_struct.Natm + 1, sizeof(int));

    // Beads are directly bonded to each other
    for (size_t i = 0; i < mol_cg.res.size(); i++)
        for (size_t j = 0; j < 3; j++) {
            _Bead bead = mol_cg.res[i].bead[j];
            if (bead.name.empty ()) continue;
            if (bead.connectedTo.size() > 0)
                for (size_t k = 0; k < bead.connectedTo.size(); k++)
                    topol_struct.con_matrix [bead.index][bead.connectedTo[k]] = 1;
        }

    // Beads are connected via another bead
    for (size_t i = 0; i < mol_cg.res.size() - 1; i++) {
        int RES = topol_struct.res_simu2res_pdb [i];
        bool ter = std::binary_search (topol_struct.ter_card.begin(), topol_struct.ter_card.end(), RES);

        /// P - P
        _Bead ph = mol_cg.res[i].bead[0];

        if (not ph.name.empty() && not ter) {
            _Bead ph2 = mol_cg.res[i+1].bead[0];
            if (not ph2.name.empty()) {
                topol_struct.con_matrix [ph.index][ph2.index] = 1;
                topol_struct.con_matrix [ph2.index][ph.index] = 1;
            }
        }

        /// RI - RI
        _Bead ri = mol_cg.res[i].bead[1];

        if (not ph.name.empty() && not ter) {
            _Bead ri2 = mol_cg.res[i+1].bead[1];
            if (not ri2.name.empty()) {
                topol_struct.con_matrix [ri.index][ri2.index] = 1;
                topol_struct.con_matrix [ri2.index][ri.index] = 1;
            }
        }

        // PH - Base
        _Bead base = mol_cg.res[i].bead[2];

        if (not base.name.empty()) {
            if (not ph.name.empty()) {
                topol_struct.con_matrix [base.index][ph.index] = 1;
                topol_struct.con_matrix [ph.index][base.index] = 1;
            }

            if (not ter) {
                _Bead ph2 = mol_cg.res[i+1].bead[0];
                if (not ph2.name.empty()) {
                    topol_struct.con_matrix [base.index][ph2.index] = 1;
                    topol_struct.con_matrix [ph2.index][base.index] = 1;
                }
            }
        }
    }

    /*for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        printf ("%4d ", i);
        for (int j = 1; j <= topol_struct.Natm_biomol; j++)
            if (topol_struct.con_matrix [i][j]) printf ("%4d ", j);
        printf ("\n");
    }*/
}

#endif
