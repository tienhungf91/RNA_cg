#include <cmath>
#include <sstream>
#include <getopt.h>
#include <dirent.h>

#include "pdb.h"
#include "mol_struct.h"

///////////////////////////////////////////////////////////////////
////////             MAIN PROGRAM                   //////////////
////////////////////////////////////////////////////////////////////
int main (int argc, char * argv []) {

    std::vector<std::vector <double>> distance;
    distance.resize(3);

    /// Read all pdb files in the current directory
    DIR *dir = NULL;
    struct dirent *file = NULL;
    std::string directory = ".";
    if ((dir = opendir (directory.c_str())) != NULL)
        while ((file = readdir (dir)) != NULL) {
            bool check = 0;
            std::string filename = file->d_name;
            size_t found = filename.find_last_of(".");
            if (filename.substr (found + 1) == "pdb") {
                printf ("Working with %s ...\n", filename.c_str());

                _Mol_aa mol_aa;
                std::vector<int> ter_card;
                Read_PDB_ion (filename, mol_aa, ter_card);
            
                _Mol_cg mol_cg = TIS_Coarse_graining_ion (mol_aa, ter_card, -1);
                //std::string pdb_out = "test.pdb";
                //Write_pdb_cg (mol_cg, pdb_out);

                for (size_t i = 0; i < mol_cg.res.size(); i++) {
                    _Res_cg res1 = mol_cg.res[i];
                    if (res1.name != "MG") continue;
                    for (size_t j = 0; j < mol_cg.res.size(); j++) {
                        _Res_cg res2 = mol_cg.res[j];
                        if (res2.name == "MG" or i == j) continue;
                        for (size_t k = 0; k < res2.bead.size(); k++) {
                            if (res2.bead[k].name.empty()) continue;
                            double dx = res1.bead[0].coord.x - res2.bead[k].coord.x;
                            double dy = res1.bead[0].coord.y - res2.bead[k].coord.y;
                            double dz = res1.bead[0].coord.z - res2.bead[k].coord.z;
                            double d = dx*dx + dy*dy + dz*dz;
                            if (d <= 100.) distance[k].push_back (sqrt(d));
                            if (d <= 16. and (k == 1 or k == 2)) check = 1;
                        }
                    }
                }
            }
            if (check) printf ("Check file %s\n", filename.c_str());
        }

    std::string file_out = "dist_Mg-P.out";
    std::ofstream OUT1 (file_out.c_str());
    if (OUT1.is_open()) {
        for (size_t i = 0; i < distance[0].size(); i++)
            OUT1 << distance[0][i] << std::endl;
        OUT1.close();
    } else {
        std::cerr << "Unable to write to file  " << file_out << std::endl;
        exit (0);
    }

    file_out = "dist_Mg-S.out";
    std::ofstream OUT2 (file_out.c_str());
    if (OUT2.is_open()) {
        for (size_t i = 0; i < distance[1].size(); i++)
            OUT2 << distance[1][i] << std::endl;
        OUT2.close();
    } else {
        std::cerr << "Unable to write to file  " << file_out << std::endl;
        exit (0);
    }

    file_out = "dist_Mg-B.out";
    std::ofstream OUT3 (file_out.c_str());
    if (OUT3.is_open()) {
        for (size_t i = 0; i < distance[2].size(); i++)
            OUT3 << distance[2][i] << std::endl;
        OUT3.close();
    } else {
        std::cerr << "Unable to write to file  " << file_out << std::endl;
        exit (0);
    }

    return (0);
}
