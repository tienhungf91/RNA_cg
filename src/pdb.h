#ifndef PDB_H
#define PDB_H

#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "DS.h"

void Assign_mass (_Atom &atom) {
    if (atom.element == "H")
        atom.mass = 1.;
    else if (atom.element == "C")
        atom.mass = 12.;
    else if (atom.element == "O")
        atom.mass = 16.;
    else if (atom.element == "N")
        atom.mass = 14.;
    else if (atom.element == "S")
        atom.mass = 32.;
    else if (atom.element == "P")
        atom.mass = 31.;
    else if (atom.element == "MG")
        atom.mass = 24.;
    else {
        std::cerr << "Unable to assign mass to atom     " << atom.name << std::endl;
        std::cerr << "             element assigned     " << atom.element << std::endl;
        exit (0);
    }
}

///////////////////////////////////////////////////////
std::string Res_type (const _Res_aa &res) {
    if ((res.name == "G") or (res.name == "C") or \
        (res.name == "U") or (res.name == "A"))
        return "RNA";
    else if ((res.name == "DG") or (res.name == "DC") or \
             (res.name == "DT") or (res.name == "DA"))
        return "DNA";
    else if ((res.name == "ALA") or (res.name == "ARG") or \
             (res.name == "ASN") or (res.name == "ASP") or \
             (res.name == "CYS") or (res.name == "GLU") or \
             (res.name == "GLN") or (res.name == "GLY") or \
             (res.name == "HIS") or (res.name == "HID") or \
             (res.name == "ILE") or (res.name == "LEU") or \
             (res.name == "LYS") or (res.name == "MET") or \
             (res.name == "PHE") or (res.name == "PRO") or \
             (res.name == "SER") or (res.name == "THR") or \
             (res.name == "TRP") or (res.name == "TYR") or \
             (res.name == "VAL"))
        return "Protein";
    else if (res.name == "MG")
        return "Ion";
    else {
        std::cerr << "Unable to recognize residue    " << res.name << std::endl;
        exit (0);
    }
}

/////////////////////////////////////////////////
void Read_PDB (const std::string &pdb_file,
               _Mol_aa &mol_aa,
               std::vector<int> &ter_card) {

    std::ifstream PDBFILE (pdb_file.c_str());
    if (PDBFILE.is_open()) {
        std::string line;
        int resid_check = -1;
        int atom_index = 0;
        _Res_aa residue, empty_residue;
        while (getline(PDBFILE, line))
            if (line.find("ATOM") == 0) {
                atom_index++;
                _Atom atom;
                atom.name = line.substr(12,4).c_str();
                // Remove spaces
                atom.name.erase (remove_if (atom.name.begin(), atom.name.end(), isspace), atom.name.end());
                atom.element = atom.name.substr(0,1).c_str();
                atom.index = atom_index;
                Assign_mass (atom);
                atom.coord.x = atof (line.substr(30,8).c_str());
                atom.coord.y = atof (line.substr(38,8).c_str());
                atom.coord.z = atof (line.substr(46,8).c_str());

                int resid_current = atoi (line.substr(22,4).c_str());
                if (resid_check == -1) resid_check = resid_current;

                if (resid_current != resid_check) {
                    // TODO: At this moment, only link consecutive residues
                    // Push the "old" residue
                    residue.connectedTo.push_back (resid_current);
                    mol_aa.res.push_back (residue);

                    // Reset everything and link the new residue to the old one
                    residue = empty_residue;
                    residue.connectedTo.push_back (resid_check);
                    resid_check = resid_current;
                }

                residue.atom.push_back (atom);
                if (residue.name.empty()) {
                    residue.index = resid_current;
                    residue.name = line.substr(17,3).c_str();
                    residue.name.erase (remove_if (residue.name.begin(), residue.name.end(), isspace), residue.name.end());
                    std::string mol_type = Res_type (residue);
                    if (mol_aa.type.empty())
                        mol_aa.type = mol_type;
                    else if (mol_aa.type != mol_type) {
                        std::cerr << "Currently only support one type of molecules in PDB file !!!!" << std::endl;
                        std::cerr << "Encounter two types   " << mol_aa.type << "     and      " << mol_type << std::endl;
                        exit (0);
                    }
                }
            } else if (line.find("TER") == 0) {
                mol_aa.res.push_back (residue);
                ter_card.push_back (residue.index);
                residue = empty_residue;
                resid_check = -1;
            }

        // Push the last residue if necessary
        if (residue.atom.size() > 0) mol_aa.res.push_back (residue);

        PDBFILE.close();
        // Sort TER card for binary search later
        std::sort (ter_card.begin(), ter_card.end());
        printf ("TER card \n");
        for (int i = 0; i < ter_card.size(); i++)
            printf ("%d\n", ter_card[i]);
    } else {
        std::cerr << "Unable to open file     " << pdb_file << std::endl;
        exit (0);
    }
}

///////////////////////////////////////////////////////////
void Read_PDB_ion (const std::string &pdb_file,
                   _Mol_aa &mol_aa,
                   std::vector<int> &ter_card) {

    std::ifstream PDBFILE (pdb_file.c_str());
    if (PDBFILE.is_open()) {
        std::string line;
        int resid_check = -1;
        int atom_index = 0;
        _Res_aa residue, empty_residue;
        while (getline(PDBFILE, line))
            if (line.find("ATOM") == 0 or line.find("HETATM") == 0) {
                atom_index++;
                _Atom atom;
                atom.name = line.substr(12,4).c_str();
                // Remove spaces
                atom.name.erase (remove_if (atom.name.begin(), atom.name.end(), isspace), atom.name.end());
                //printf ("%s\n", atom.name.c_str());
                atom.element = atom.name.substr(0,1).c_str();

                atom.index = atom_index;
                atom.coord.x = atof (line.substr(30,8).c_str());
                atom.coord.y = atof (line.substr(38,8).c_str());
                atom.coord.z = atof (line.substr(46,8).c_str());

                int resid_current = atoi (line.substr(22,4).c_str());
                if (resid_check == -1) resid_check = resid_current;

                if (resid_current != resid_check) {
                    // TODO: At this moment, only link consecutive residues
                    // Push the "old" residue
                    residue.connectedTo.push_back (resid_current);
                    mol_aa.res.push_back (residue);

                    // Reset everything and link the new residue to the old one
                    residue = empty_residue;
                    residue.connectedTo.push_back (resid_check);
                    resid_check = resid_current;
                }

                //Assign_mass (atom);
                residue.atom.push_back (atom);

                if (residue.name.empty()) {
                    residue.index = resid_current;
                    residue.name = line.substr(17,3).c_str();
                    residue.name.erase (remove_if (residue.name.begin(), residue.name.end(), isspace), residue.name.end());
                    std::string mol_type = Res_type (residue);
                    if (mol_aa.type.empty() and mol_type != "Ion")
                        mol_aa.type = mol_type;
                    else if (mol_aa.type != mol_type and mol_type != "Ion" and mol_aa.type != "Ion") {
                        std::cerr << "Currently only support one type of molecules in PDB file !!!!" << std::endl;
                        std::cerr << "Encounter two types   " << mol_aa.type << "     and      " << mol_type << std::endl;
                        exit (0);
                    }
                }
            } else if (line.find("TER") == 0) {
                mol_aa.res.push_back (residue);
                ter_card.push_back (residue.index);
                residue = empty_residue;
                resid_check = -1;
            }

        // Push the last residue if necessary
        if (residue.atom.size() > 0) mol_aa.res.push_back (residue);

        PDBFILE.close();
        // Sort TER card for binary search later
        std::sort (ter_card.begin(), ter_card.end());
        printf ("TER card \n");
        for (int i = 0; i < ter_card.size(); i++)
            printf ("%d\n", ter_card[i]);
    } else {
        std::cerr << "Unable to open file     " << pdb_file << std::endl;
        exit (0);
    }
}

/////////////////////////////////////////////////
void provide_coord (const std::string &pdb_file,
                    const _topol_struct &topol_struct,
                    _coord_vel_force_struct &coord_vel_force_struct) {

    std::ifstream PDBFILE (pdb_file.c_str());
    if (PDBFILE.is_open()) {
        std::string line;
        while (getline(PDBFILE, line))
            if (line.find("ATOM") == 0) {
                int RES = atoi (line.substr(22,4).c_str());
                int res = topol_struct.res_pdb2res_simu [RES];

                std::string bead_name = line.substr (12,4).c_str();
                bead_name.erase (remove_if (bead_name.begin(), bead_name.end(), isspace), bead_name.end());
                int atm_key = 0;
                if (bead_name == "PH")
                    atm_key = topol_struct.atom_key[res][0];
                else if (bead_name == "RI")
                    atm_key = topol_struct.atom_key[res][1];
                else atm_key = topol_struct.atom_key[res][2];

                /// Parse coord
                coord_vel_force_struct.coordx [atm_key] = atof (line.substr(30,8).c_str());
                coord_vel_force_struct.coordy [atm_key] = atof (line.substr(38,8).c_str());
                coord_vel_force_struct.coordz [atm_key] = atof (line.substr(46,8).c_str());
            }

        PDBFILE.close();
    } else {
        std::cerr << "Unable to open file     " << pdb_file << std::endl;
        exit (0);
    }
}

////////////////////////////////////////////////////
void Read_PDB_cg (const std::string &pdb_file,
                  int &Natm_biomol, int &Natm,
                  _Mol_cg &mol_cg) {

    Natm_biomol = 0;
    Natm   = 0;

    std::ifstream PDBFILE (pdb_file.c_str());
    if (PDBFILE.is_open()) {
        std::string line;
        int resid_check = -1;
        _Res_cg residue, empty_residue;
        while (getline(PDBFILE, line))
            if (line.find("ATOM") == 0) {
                Natm++;
                _Bead bead;
                bead.name = line.substr(12,4).c_str();
                // Remove spaces
                bead.name.erase (remove_if (bead.name.begin(), bead.name.end(), isspace), bead.name.end());
                if (bead.name == "Mg") {
               //     N_ion1++;
                    continue;
                } else if (bead.name == "Cl") {
               //     N_ion2++;
                    continue;
                } else if (bead.name == "K") {
               //     N_ion3++;
                    continue;
                } else Natm_biomol++;

                bead.index = atoi (line.substr(6,5).c_str());
             //   bead.coord.x = atof (line.substr(30,8).c_str());
             //   bead.coord.y = atof (line.substr(38,8).c_str());
             //   bead.coord.z = atof (line.substr(46,8).c_str());

                size_t resid_current = atoi (line.substr(22,4).c_str());
                if (resid_check == 0) resid_check = resid_current;

                if (resid_current != resid_check) {
                    // TODO: At this moment, only link consecutive residues
                    // Push the "old" residue
                    residue.connectedTo.push_back (resid_current);
                    mol_cg.res.push_back (residue);

                    // Reset everything and link the new residue to the old one
                    residue = empty_residue;
                    residue.connectedTo.push_back (resid_check);
                    resid_check = resid_current;
                }

                residue.bead.push_back (bead);
                if (residue.name.empty()) {
                    residue.index = resid_current;
                    residue.name = line.substr(17,3).c_str();
                    residue.name.erase (remove_if (residue.name.begin(), residue.name.end(), isspace), residue.name.end());
                  /*  std::string mol_type = Res_type (residue);
                    if (mol_aa.type.empty())
                        mol_aa.type = mol_type;
                    else if (mol_aa.type != mol_type) {
                        std::cerr << "Currently only support one type of molecules in PDB file !!!!" << std::endl;
                        std::cerr << "Encounter two types   " << mol_aa.type << "     and      " << mol_type << std::endl;
                        exit (0);
                    }*/
                }
            }/* else if (line.find("TER") == 0) {
                mol_aa.res.push_back (residue);
                ter_card.push_back (residue.index);
                residue = empty_residue;
                resid_check = 0;
            }*/

        // Push the last residue if necessary
        if (residue.bead.size() > 0) mol_cg.res.push_back (residue);

        PDBFILE.close();
    } else {
        std::cerr << "Unable to open file     " << pdb_file << std::endl;
        exit (0);
    }
}

//////////////////////////////////////////////////////////////////////
void Write_pdb_cg (const _Mol_cg &mol_cg,
                   const std::string &pdb_out) {
    std::ofstream OUT (pdb_out.c_str());
    if (OUT.is_open()) {
        for (size_t i = 0; i < mol_cg.res.size(); i++)
            for (size_t j = 0; j < mol_cg.res[i].bead.size(); j++) {
                _Bead bead = mol_cg.res[i].bead[j];
                if (bead.name.empty ()) continue;
                std::string beadname;
                if (bead.name.length() < 4) {
                    beadname = " " + bead.name;
                    for (size_t j = 0; j < 3 - bead.name.length(); j++)
                        beadname = beadname + " ";
                }
                else beadname = bead.name;
                OUT << "ATOM  " << std::setw(5) << bead.index << " " << std::setw(4) << std::left << beadname;
                OUT << std::setw(4) << std::right << mol_cg.res[i].name << std::setw(6) << mol_cg.res[i].index \
                    << std::setw(12) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << bead.coord.x \
                    << std::setw(8) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << bead.coord.y \
                    << std::setw(8) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << bead.coord.z << std::endl;
            }

        /// Write out connectivity
        OUT << "TER\n";
        for (size_t i = 0; i < mol_cg.res.size(); i++)
            for (size_t j = 0; j < mol_cg.res[i].bead.size(); j++) {
                _Bead bead = mol_cg.res[i].bead[j];
                if (bead.name.empty ()) continue;
                if (bead.connectedTo.size() > 0) {
                    OUT << "CONECT" << std::setw(5) << bead.index;
                    for (size_t j = 0; j < bead.connectedTo.size(); j++)
                        OUT << std::setw (5) << bead.connectedTo [j];
                    OUT << std::endl;
                }
            }
        OUT.close();
    }
    else {
        std::cerr << "Unable to write to file   " << pdb_out << std::endl;
        exit (0);
    }
}

////////////////////////////////////////////////////////////
void Write_pdb_cg2 (const _Mol_cg &mol_cg,
                    const _topol_struct &topol_struct,
                    const _crowder_struct &crowder_struct,
                    const _coord_vel_force_struct &coord_vel_force_struct,
                    int ion_type,
                    const std::string &pdb_out) {
    std::ofstream OUT (pdb_out.c_str());
    if (OUT.is_open()) {
        int index = 0;
        for (size_t i = 0; i < mol_cg.res.size(); i++)
            for (size_t j = 0; j < 3; j++) {
                _Bead bead = mol_cg.res[i].bead[j];
                if (bead.name.empty ()) continue;
                index++;
                std::string beadname;
                if (bead.name.length() < 4) {
                    beadname = " " + bead.name;
                    for (size_t j = 0; j < 3 - bead.name.length(); j++)
                        beadname = beadname + " ";
                }
                else beadname = bead.name;
                OUT << "ATOM  " << std::setw(5) << bead.index << " " << std::setw(4) << std::left << beadname;
                OUT << std::setw(4) << std::right << mol_cg.res[i].name << std::setw(6) << mol_cg.res[i].index \
                    << std::setw(12) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << coord_vel_force_struct.coordx [index] \
                    << std::setw(8) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << coord_vel_force_struct.coordy [index] \
                    << std::setw(8) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << coord_vel_force_struct.coordz [index] << std::endl;
            }

        for (int i = topol_struct.Natm_biomol + 1; i <= topol_struct.Natm; i++) {
            std::string name;
            if (i <= topol_struct.Natm_biomol + crowder_struct.N_crwd[0])
                if (ion_type == 1)
                    name = "Mg";
                else if (ion_type == 2)
                    name = "Ca";
                else if (ion_type == 3)
                    name = "Co";
                else {
                    printf ("Not recognize ion type!!!!\n");
                    exit (0);
                }
            else if (i <= topol_struct.Natm_biomol + crowder_struct.N_crwd[0] + crowder_struct.N_crwd[1])
                name = "Cl";
            else
                name = "K ";

            OUT << "ATOM  " << std::setw(5) << i << "  " << std::setw(4) << std::left << name;
            OUT << std::setw(4) << std::right << name << std::setw(5) << i \
                << std::setw(12) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << coord_vel_force_struct.coordx [i] \
                << std::setw(8) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << coord_vel_force_struct.coordy [i] \
                << std::setw(8) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << coord_vel_force_struct.coordz [i] << std::endl;
        }

        /// Write out connectivity
        OUT << "TER\n";
        for (size_t i = 0; i < mol_cg.res.size(); i++)
            for (size_t j = 0; j < 3; j++) {
                _Bead bead = mol_cg.res[i].bead[j];
                if (bead.name.empty ()) continue;
                if (bead.connectedTo.size() > 0) {
                    OUT << "CONECT" << std::setw(5) << bead.index;
                    for (size_t k = 0; k < bead.connectedTo.size(); k++)
                        OUT << std::setw (5) << bead.connectedTo [k];
                    OUT << std::endl;
                }
            }
        OUT.close();
    }
    else {
        std::cerr << "Unable to write to file   " << pdb_out << std::endl;
        exit (0);
    }
}

////////////////////////////////////////////////////////////
void Write_pdb_aa (const _Mol_aa &mol_aa,
                   const std::string &pdb_out) {
    std::ofstream OUT (pdb_out.c_str());
    if (OUT.is_open()) {
        for (size_t i = 0; i < mol_aa.res.size(); i++) {
            _Res_aa res = mol_aa.res[i];
            for (size_t j = 0; j < res.atom.size(); j++) {
                _Atom atom = res.atom[j];
                std::string atomname;
                if (atom.name.length() < 4) {
                    atomname = " " + atom.name;
                    for (size_t j = 0; j < 3 - atom.name.length(); j++)
                        atomname = atomname + " ";
                } else atomname = atom.name;

                OUT << "ATOM  " << std::setw(5) << atom.index << " " << std::setw(4) << std::left << atomname;
                OUT << std::setw(4) << std::right << res.name << std::setw(6) << res.index \
                    << std::setw(12) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << atom.coord.x \
                    << std::setw(8) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << atom.coord.y \
                    << std::setw(8) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << atom.coord.z << std::endl;
            }
        }
        OUT.close();
    } else {
        std::cerr << "Unable to write to file   " << pdb_out << std::endl;
        exit (0);
    }
}

#endif
