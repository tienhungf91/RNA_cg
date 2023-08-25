/* Unit used in the code
   Mass     m0 = 1 Da = 1.6605 x 10^(-27) kg
   Distance a0 = 1 Angstrom
   Energy   e0 = 1 kcal/mol = 6.955 x 10^(-21) kg.m2/s2
   Time     t0 = sqrt (m0/e0) * a0 = 4.886336 x 10^(-14) s ~ 50 fs
   Temperature T0 = kT/e0
   Viscosity of water n(water) ~ 10^(-3) Pa.s
   In low friction regime:
       n_L = 0.01 n(water) ~ 0.029427 Da.A^(-1) t0^(-1)
       gamma = 6PI * n_L * Ri ~ 0.555 Ri Da.A^(-1) t0^(-1)
       Time step = dt * t0 = 0.05 * 50 fs
   In high friction regime:
       gamma ~ 55.5 Ri
       Time step = 0.5 * 50 fs
   All the time steps above are just examples, need to perform trial-and-error run to determine appropriate values
*/

#include "dcd.h"
#include <cmath>
#include <sstream>
#include <getopt.h>

#include "simpleMath.h"
#include "pdb.h"
#include "mol_struct.h"
#include "crowder.h"
#include "force.h"
#include "hydbond.h"
#include "stack.h"
#include "neighborlist.h"
#include "mani_coord.h"
#include "IO.h"
#include "debug.h"
#include "mem.h"
#include "EW.h"
#include "analysis.h"
#include "time.h"

////////////////////////////////////////////////////////////////////////
void usage () {
    printf("        Running options:\n");
    printf("   -p     --pdb                pdb file\n");
    printf("   -n     --pdb_CG             if provided, use coordinates from this pdb file\n");
    printf("   -l     --ligand             ligand index\n");
    printf("   -P     --out_pdb            output pdb file for visualization   [ system.pdb ]\n");
    printf("   -b     --bond_file          hydrogen bond file\n");
    printf("   -k     --stack_file         tertiary stack file\n");
    printf("   -m     --maxi_file          maxi key file\n");
    printf("   -t     --traj_file          trajectory file    [ md.dcd ]\n");
    printf("   -o     --mdout              md output file     [ md.out ]\n");
    printf("   -i     --info_file          info file          [ mdinfo ]\n");
    printf("   -r     --rst_file           restart file       [ md.rst ]\n");
    printf("   -u     --uvv_file           ion-P pair potential file\n");
    printf("   -T     --temp               temperature (oC)\n");
    printf("   -D     --ion_type           specify which ion used (1=Mg [default]; 2=Ca; 3=CoHex)\n");
    printf("   -M     --Mgconc             di or trivalent ion concentration\n");
    printf("   -K     --Kconc              monovalent ion concentration\n");
    printf("   -s     --nstep              number of step\n");
    printf("   -h     --timestep           timestep [default 0.05]\n");
    printf("   -g     --gamma              viscosity [default 0.555]\n");
    printf("   -v     --box                box size (x y z)\n");
    printf("   -a     --step_traj          step to write restart and trajectory files\n");
    printf("   -e     --step_ener          step to write energy to output file\n");
//    printf("   -B     --length_per_charge  length per unit (bare) charge in the RNA [4.4 A]\n");
//    printf("   -G     --deltaG0            free energy correction of 2nd stack [0.80 kcal/mol]\n");
//    printf("   -U     --Ust0               3rd stack energy [6.5 kcal/mol]\n");
//    printf("   -H     --Uhb0               hydrogen bond energy [2.70 kcal/mol]\n");
//    printf("   -w     --weight_NN          weight for non-native interactions [1.0]\n");
    printf("   -c     --fix_solute         flag to fix the solute\n");
    printf("   -C     --Cutoff             cutoff distance for electrostatic force [20 A]\n");
    printf("   -R     --restart            flag to restart\n");
    printf("   -E     --Equilibrate        number of step to equilibrate ions first\n");
//    printf("   -d     --buffer             buffer distance (for neighbor-list update) [3 A]\n");
    printf("   -f     --force              constant force (pN) in pulling simulation\n");
    printf("   -F     --f_direction        force direction [default \"0 0 1\"]\n");
    printf("   -I     --f_id               two bead IDs to apply the external force (id1 id2)\n");
    printf("   -z     --zeroQ              set those phosphate charges to 0\n");
    printf("   -S     --seed               random seed\n");
//#ifdef _OPENMP
//    printf("   -n     --ncpu               number of cpus used\n");
//#endif
}

////////////////////////////////////////
/// Calculate temperature (for testing)
////////////////////////////////////////
/*double calc_temp (const _topol_struct           &topol_struct,
                  const _coord_vel_force_struct &coord_vel_force_struct) {
    double temp = 0;
    for (int i = 1; i <= topol_struct.Natm; i++) {
        double x = coord_vel_force_struct.velox [i] * coord_vel_force_struct.velox [i];
        double y = coord_vel_force_struct.veloy [i] * coord_vel_force_struct.veloy [i];
        double z = coord_vel_force_struct.veloz [i] * coord_vel_force_struct.veloz [i];
        int j = topol_struct.maxi_key [i];
        temp += (x + y + z) * topol_struct.MASS [j];
    }
    return ONE_THIRD * temp / (topol_struct.Natm * KELVIN_TO_KT);
}*/

///////////////////////////////////////////////////////////////
void deterministic_forces (const std::vector<_Res_cg> &res,
                           _simu_struct            &simu_struct,
                           _topol_struct           &topol_struct,
                           _hydbond_struct         &hydbond_struct,
                           _stack_struct           &stack_struct,
                           _coord_vel_force_struct &coord_vel_force_struct,
                           _energy_struct          &energy_struct,
                           unsigned int            *myseed) {
//                           const _EW_struct        &EW_struct) {

    energy_struct.E_BOND = 0;
    energy_struct.E_ANGLE = 0;
    energy_struct.E_EXCLUDED = 0;
    energy_struct.E2_HB = 0;
    energy_struct.E3_HB = 0;
    energy_struct.E2_STACK = 0;
    energy_struct.E3_STACK = 0;
    energy_struct.E_Q = 0;

/*    if (simu_struct.relist_step % simu_struct.LR_skip == 0) {
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
    SR_interactions_lists (simu_struct, topol_struct, coord_vel_force_struct, EW_struct);*/

    for (int i = 1; i <= topol_struct.Natm; i++) {
        coord_vel_force_struct.forcex [i] = 0;
        coord_vel_force_struct.forcey [i] = 0;
        coord_vel_force_struct.forcez [i] = 0;
    }
    non_bonded_interaction (simu_struct, topol_struct, coord_vel_force_struct, energy_struct);
    polymer_constraints (simu_struct.box, res, topol_struct, coord_vel_force_struct, energy_struct);
    stacking_interactions (topol_struct, coord_vel_force_struct, stack_struct, energy_struct);
    hydrogen_bonds_interactions (coord_vel_force_struct, hydbond_struct, simu_struct.temp, energy_struct, myseed);
}



///////////////////////////////////////////////////////////////////
////////             MAIN PROGRAM                   //////////////
////////////////////////////////////////////////////////////////////
int main (int argc, char * argv []) {

    _simu_struct             simu_struct;
    _topol_struct            topol_struct;
    _coord_vel_force_struct  coord_vel_force_struct;
    _crowder_struct          crowder_struct;
    _stack_struct            stack_struct;
    _hydbond_struct          hydbond_struct;
    _energy_struct           energy_struct;
//    _EW_struct               EW_struct;

//    EW_struct.sigma          = 20.0;   // sigma < box length / 5
//    EW_struct.kappa          = 1. / EW_struct.sigma;
//    EW_struct.CTM            = 2500;
//    EW_struct.N_vec          = 0;
//    EW_struct.EXP_CUTOFF     = 5.0;
//    EW_struct.sinus          = NULL;
//    EW_struct.kosinus        = NULL;

    simu_struct.step_traj    = 50000;
    simu_struct.step_energy  = 10000;
    simu_struct.nstep        = 100000000;
    simu_struct.restraint_step = 0;   // Perform MD with the biomolecule fixed
    simu_struct.step_adj_box = 100;
    simu_struct.visc_gamma   = 0.555;
    simu_struct.dt           = 0.05;
//    simu_struct.LR_skip      = 1;     // Skip LR interactions for every LR_skip step
    simu_struct.progress     = 0;
    simu_struct.restart      = 0;
    simu_struct.fix_solute   = 0;
    simu_struct.debug        = 0;

    simu_struct.pull_force   = 0;
    simu_struct.f_direction.x = 0;    simu_struct.f_direction.y = 0;     simu_struct.f_direction.z = 1;
    simu_struct.force        = 0;

    simu_struct.Coulomb      = 0;
//    simu_struct.Ze           = 1;
//    simu_struct.epsilon      = 78;
    simu_struct.b            = 4.38178046;  // length per unit (bare) charge
    simu_struct.kappaD       = 0;
//    simu_struct.lambdaD      = 9e9;
//    simu_struct.D_SR_CUTOFF  = 12.0;
    simu_struct.D_CUTOFF  = 20.0;
//    simu_struct.D_CT_CUTOFF  = 1.0 * EW_struct.sigma;
//    simu_struct.D2_CT_CUTOFF = simu_struct.D_CT_CUTOFF * simu_struct.D_CT_CUTOFF;
//    simu_struct.D2_CUTOFF = simu_struct.D_CUTOFF * simu_struct.D_CUTOFF;
    simu_struct.Mc_total     = 0;
    simu_struct.cell_content = NULL;
    simu_struct.relist       = 0.;
//    simu_struct.relist_step  = 0;
    simu_struct.dr1          = 3.;     // buffer distance for neighbor list updating
    simu_struct.ion_type     = 1;      // Mg is the default
    simu_struct.ion_q        = 2;

   // coord_vel_force_struct.coordx     = NULL;
   // coord_vel_force_struct.coordy     = NULL;
   // coord_vel_force_struct.coordz     = NULL;
   // coord_vel_force_struct.coord_oldx = NULL;
   // coord_vel_force_struct.coord_oldy = NULL;
   // coord_vel_force_struct.coord_oldz = NULL;
   // coord_vel_force_struct.velox      = NULL;
   // coord_vel_force_struct.veloy      = NULL;
   // coord_vel_force_struct.veloz      = NULL;
    coord_vel_force_struct.forcex     = NULL;
    coord_vel_force_struct.forcey     = NULL;
    coord_vel_force_struct.forcez     = NULL;
//    coord_vel_force_struct.flx        = NULL;
//    coord_vel_force_struct.fly        = NULL;
//    coord_vel_force_struct.flz        = NULL;

    topol_struct.N_amino     = 0;
    topol_struct.Nres_biomol = 0;
    topol_struct.Nres        = 0;
    topol_struct.Natm_biomol = 0;
    topol_struct.Npair_biomol = 0;
    topol_struct.Natm        = 0;
    topol_struct.bead_center = -1;
    topol_struct.atom_key    = NULL;
    topol_struct.ligand_index = -1;
    topol_struct.ligand_bead = -1;

    stack_struct.st_D        = 5.0;   // Tertiary stack U_st
    stack_struct.ss_D        = 0.9;  // Secondary stack deltaG_0
    stack_struct.s3_N        = 0;
    stack_struct.st_N        = 0;
    stack_struct.score_cutoff = 0.2;
//    stack_struct.ST_ATOM     = NULL;

    hydbond_struct.hs_D      = 2.7;  // Hbond U_HB
    hydbond_struct.hb_N      = 0;
    hydbond_struct.hb_code   = NULL;
    hydbond_struct.HB_Natm   = 0;
    hydbond_struct.HB_atom_key = NULL;
    hydbond_struct.HB_ATOM   = NULL;
    hydbond_struct.ATOM_HB   = NULL;
    hydbond_struct.HB_PAIR   = NULL;
    hydbond_struct.weight_nonnative = 1.0;

    energy_struct.E_BOND     = 0;
    energy_struct.E_ANGLE    = 0;
    energy_struct.E_EXCLUDED = 0;
    energy_struct.E2_HB      = 0;
    energy_struct.E3_HB      = 0;
    energy_struct.E2_STACK   = 0;
    energy_struct.E3_STACK   = 0;
    energy_struct.E_Q        = 0;

    char * file_PDB = NULL;
    char * file_PDB_sup = NULL;
    char * file_unprocessed_bonds  = NULL;
    char * file_unprocessed_stacks = NULL;
    char * file_maxi;
    char * file_dcd        = "md.dcd";      //// Trajectory file
    char * file_output     = "md.out";      //// Output energy and such
    char * file_info       = "mdinfo";      //// Reading progress and restart parameter
    char * file_rst        = "md.rst";      //// Similar format as AMBER rst file
                                            //// Contain coords and velocities for restart
    std::string out_pdb    = "system.pdb";  //// pdb output, for visualizing purpose
    std::string uvv_file;
    std::string temp_string, zeroQ;

//    char * file_error      = "Error.dat";

    int option_char;
    do {
        struct option long_options[] =
        {
            {"pdb",                required_argument,  0,  'p'},
            {"pdb_CG",             required_argument,  0,  'n'},
            {"ligand",             required_argument,  0,  'l'},
            {"out_pdb",            required_argument,  0,  'P'},
            {"bond_file",          required_argument,  0,  'b'},
            {"stack_file",         required_argument,  0,  'k'},
            {"maxi_key",           required_argument,  0,  'm'},
            {"traj_file",          required_argument,  0,  't'},
            {"mdout",              required_argument,  0,  'o'},
            {"info_file",          required_argument,  0,  'i'},
            {"rst_file",           required_argument,  0,  'r'},
            {"uvv_file",           required_argument,  0,  'u'},
            {"temp",               required_argument,  0,  'T'},
            {"ion_type",           required_argument,  0,  'D'},
            {"Mgconc",             required_argument,  0,  'M'},
            {"Kconc",              required_argument,  0,  'K'},
            {"box",                required_argument,  0,  'v'},
            {"nstep",              required_argument,  0,  's'},
            {"timestep",           required_argument,  0,  'h'},
            {"gamma",              required_argument,  0,  'g'},
            {"fix_solute",         no_argument,        0,  'c'},
//            {"length_per_charge",  required_argument,  0,  'B'},
            {"deltaG0",            required_argument,  0,  'G'},
            {"Ust0",               required_argument,  0,  'U'},
            {"Uhb0",               required_argument,  0,  'H'},
            {"weight_NN",          required_argument,  0,  'w'},
            {"Cutoff",             required_argument,  0,  'C'},
            {"restart",            no_argument,        0,  'R'},
            {"Equilibrate",        required_argument,  0,  'E'},
            {"force",              required_argument,  0,  'f'},
            {"f_direction",        required_argument,  0,  'F'},
            {"f_id",               required_argument,  0,  'I'},
            {"buffer",             required_argument,  0,  'd'},
            {"step_traj",          required_argument,  0,  'a'},
            {"step_ener",          required_argument,  0,  'e'},
            {"seed",               required_argument,  0,  'S'},
            {"zeroQ",              required_argument,  0,  'z'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_char = getopt_long (argc, argv, "p:n:l:P:b:k:m:t:o:i:r:u:T:D:M:K:v:s:h:g:cG:U:H:w:C:RE:f:F:I:d:a:e:S:z:", long_options, &option_index);
        // followed by 1 colon - require an argument; no colon - no argument required

        if (option_char == -1)
            usage();
        else switch (option_char) {
            char *end1, *end2, *end3;
            case '?':
                usage();
                exit (0);
            case 'p':
                file_PDB = optarg;
                break;
            case 'n':
                file_PDB_sup = optarg;
                break;
            case 'l':
                topol_struct.ligand_index = atoi (optarg);
                break;
            case 'P':
                out_pdb = optarg;
                break;
            case 'b':
                file_unprocessed_bonds = optarg;
                break;
            case 'k':
                file_unprocessed_stacks = optarg;
                break;
            case 'm':
                file_maxi = optarg;
                break;
            case 't':
                file_dcd = optarg;
                break;
            case 'o':
                file_output = optarg;
                break;
            case 'i':
                file_info = optarg;
                break;
            case 'r':
                file_rst = optarg;
                break;
            case 'u':
                uvv_file = optarg;
                break;
            case 'T':
                simu_struct.temp = (atof (optarg) + 273.15) * KELVIN_TO_KT;
                break;
            case 'D':
                simu_struct.ion_type = atoi (optarg);
                if (simu_struct.ion_type == 3) simu_struct.ion_q = 3;
                break;
            case 'M':
                crowder_struct.conc [0] = atof (optarg);
                break;
            case 'K':
                crowder_struct.conc [2] = atof (optarg);
                break;
            case 'v':
                simu_struct.box[0] = strtod (optarg, &end1);
                simu_struct.box[1] = strtod (end1, &end2);
                simu_struct.box[2] = strtod (end2, &end3);
                break;
            case 's':
                simu_struct.nstep = atol (optarg);
                break;
            case 'h':
                simu_struct.dt = atof (optarg);
                break;
            case 'g':
                simu_struct.visc_gamma = atof (optarg);
                break;
            case 'c':
                simu_struct.fix_solute = 1;
                break;
//            case 'B':
//                simu_struct.b = atof (optarg);
//                break;
            case 'G':
                stack_struct.ss_D = atof (optarg);
                break;
            case 'U':
                stack_struct.st_D = atof (optarg);
                break;
            case 'H':
                hydbond_struct.hs_D = atof (optarg);
                break;
            case 'w':
                hydbond_struct.weight_nonnative = atof (optarg);
                break;
            case 'C':
//                simu_struct.Coulomb = 1;
                simu_struct.D_CUTOFF = atof (optarg);
                break;
            case 'R':
                simu_struct.restart = 1;
                break;
            case 'E':
                simu_struct.restraint_step = atol (optarg);
                break;
            case 'f':
                simu_struct.pull_force = 1;
                simu_struct.force = atof (optarg);
                break;
            case 'F':
                simu_struct.f_direction.x = strtod (optarg, &end1);
                simu_struct.f_direction.y = strtod (end1, &end2);
                simu_struct.f_direction.z = strtod (end2, &end3);
                break;
            case 'I':
                simu_struct.f_id1 = strtol (optarg, &end1, 10); // base 10
                simu_struct.f_id2 = strtol (end1,   &end2, 10);
                break;
            case 'd':
                simu_struct.dr1 = atof (optarg);
                break;
            case 'a':
                simu_struct.step_traj = atoi (optarg);
                break;
            case 'e':
                simu_struct.step_energy = atoi (optarg);
                break;
//            case 'd':
//                temp_string = optarg;
//                break;
            case 'S':
                simu_struct.seed = atoi (optarg);
                break;
            case 'z':
                zeroQ = optarg;
                break;
        }
    } while (option_char != -1);

/*#ifdef _OPENMP
    omp_set_dynamic(0);
    if (ncpu == 0)
        ncpu = omp_get_max_threads ();
    omp_set_num_threads (ncpu);
#endif*/

    std::stringstream stream (zeroQ);
    int tmp_zeroQ;
    std::vector<int> res_zeroQ;
    while (stream >> tmp_zeroQ)
        res_zeroQ.push_back (tmp_zeroQ);
    std::sort (res_zeroQ.begin(), res_zeroQ.end());

    unsigned myseed = simu_struct.temp * 1000 * simu_struct.seed;
    srand (myseed);

    double t2 = simu_struct.temp * simu_struct.temp;
    simu_struct.epsilon = 296.0736276 - 619.2813716 * simu_struct.temp + 531.2826741 * t2 - 180.0369914 * simu_struct.temp * t2;
//    printf ("epsilon  =  %f\n", simu_struct.epsilon);

    simu_struct.lB_T = 332.0637090 / simu_struct.epsilon;
//    simu_struct.lB_T = 332.0637090 * (3.2736e-2 * simu_struct.temp - 6.6962575e-3); //expanded at T = 50 oC, in kcal/mol

//    simu_struct.lB_T = 332.0637090 * (3.464665e-2 * simu_struct.temp - 7.75551e-3); // linear fit in T=[0..120] oC range
    printf ("lB   = %f\n", simu_struct.lB_T / simu_struct.temp);

    simu_struct.Ze = simu_struct.b * simu_struct.temp / simu_struct.lB_T;
    printf ("Ze   = %f\n", simu_struct.Ze);
    if (crowder_struct.conc [0] > 1e-6) {
        read_pair_potential (uvv_file, simu_struct);
        update_Ze (simu_struct, crowder_struct.conc);
        printf ("new Ze   = %f\n", simu_struct.Ze);
    }

    double rho = 2*crowder_struct.conc [2] * 6.022e-4;
//    printf ("KCl concentration  %f\n", crowder_struct.conc [2]);
    simu_struct.kappaD = sqrt (4*PI * simu_struct.lB_T * rho / simu_struct.temp);
    simu_struct.lambdaD = 1. / simu_struct.kappaD;

    printf ("lambda  = %f\n", simu_struct.lambdaD);
    if (crowder_struct.conc [0] > 1e-6) update_pmf_ion_P_implicit (simu_struct.pair_potential, simu_struct);

//    simu_struct.D_CUTOFF = 3. * simu_struct.lambdaD;
//    simu_struct.D_CUTOFF = 50.;
    if (simu_struct.D_CUTOFF >= 0.5 * simu_struct.box[0] || \
        simu_struct.D_CUTOFF >= 0.5 * simu_struct.box[1] || \
        simu_struct.D_CUTOFF >= 0.5 * simu_struct.box[2]) {
        printf ("Box is too small, it should be at least %6.1f A!!\n", 2*simu_struct.D_CUTOFF);
        exit (0);
    }
    printf ("Cutoff  %6.1f A\n", simu_struct.D_CUTOFF);
    simu_struct.D2_CUTOFF = simu_struct.D_CUTOFF * simu_struct.D_CUTOFF;

//    simu_struct.lB_T_D2 = simu_struct.lB / simu_struct.D2_CUTOFF;

    /////////// HARD-CODE NS
    //////////// Not seeing topol_struct.PD anywhere else in the code
/*    for (int In1 = 0; In1 < NS; In1++)
        for (int In2 = In1; In2 < NS; In2++)
            {
                size_t In12 = NS * In1 - In1 * (In1 + 1) / 2 + In2;
                topol_struct.PD [ In12 ]    = ( topol_struct.D [ In1 ] + topol_struct.D [ In2 ] ) / 2;
                topol_struct.PD_sq [ In12 ] = topol_struct.PD [ In12 ] * topol_struct.PD [ In12 ];
            }*/

    if (file_maxi != NULL)
        read_maxi_key (file_maxi, simu_struct, topol_struct);
    else {
        printf ("Need maxikey file !!!\n");
        exit (0);
    }

    _Mol_aa mol_aa;
    if (file_PDB != NULL)
        Read_PDB (file_PDB, mol_aa, topol_struct.ter_card);
    else {
        printf ("Need a pdb file !!!\n");
        exit (0);
    }
    _Mol_cg mol_cg = TIS_Coarse_graining (mol_aa, topol_struct.ter_card, topol_struct.ligand_index);

    Parse_mol_struct_to_simu (mol_cg, topol_struct, res_zeroQ, coord_vel_force_struct);
    if (topol_struct.bead_zeroQ.size() > 0) {
        printf ("Bead having zero charge   ");
        for (int i = 0; i < topol_struct.bead_zeroQ.size(); i++)
            printf ("%6d", topol_struct.bead_zeroQ[i]);
        printf ("\n");
    }

    /*printf ("!!!!!connect!!!!!\n");
    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        printf ("%4d | ", i);
        for (int j = 1; j <= topol_struct.Natm_biomol; j++)
            if (topol_struct.con_matrix [i][j] == 1)
                printf ("%4d", j);
        printf ("\n");
    }*/

//    topol_struct.INDX = (int *) calloc (topol_struct.Natm_biomol + 1, sizeof(int));
//    topol_struct.JNDX = (int *) calloc (topol_struct.Natm_biomol + 1, sizeof(int));
    if (not simu_struct.fix_solute) {
        hydbond_struct.eval = 1;
        initialize_HB_atoms (mol_cg, topol_struct, hydbond_struct);
        if (file_unprocessed_bonds != NULL)
            initialize_unprocessed_bonds (file_unprocessed_bonds, mol_cg, topol_struct, hydbond_struct, coord_vel_force_struct);

        initialize_unprocessed_stacks (file_unprocessed_stacks, mol_aa, mol_cg, topol_struct, stack_struct, coord_vel_force_struct, simu_struct.temp);
    }

    if (file_PDB_sup != NULL)
        provide_coord (file_PDB_sup, topol_struct, coord_vel_force_struct);

    simu_struct.volume = simu_struct.box[0] * simu_struct.box[1] * simu_struct.box[2];
    list_crowder_types (simu_struct.ion_type, crowder_struct);
    crowder_struct.N_crwd [0] = (int) (crowder_struct.conc[0] * 6.022e-04 * simu_struct.volume);
//    crowder_struct.N_crwd [2] = (int) (crowder_struct.conc[2] * 6.022e-04 * simu_struct.volume);
//    crowder_struct.N_crwd [1] = 2 * crowder_struct.N_crwd [0] + crowder_struct.N_crwd [2] + mol_cg.netcharge;

//    printf ("No of K     %d\n", crowder_struct.N_crwd [2]);
//    printf ("No of Cl    %d\n", crowder_struct.N_crwd [1]);
//    topol_struct.Nres = topol_struct.Nres_biomol + crowder_struct.N_crwd [0] + crowder_struct.N_crwd [1] + crowder_struct.N_crwd [2];
    topol_struct.Nres = topol_struct.Nres_biomol + crowder_struct.N_crwd [0];
    allocate_mem (topol_struct, crowder_struct, coord_vel_force_struct);
    
    if (simu_struct.restart == 0) {
        ///// Align RNA axes and move RNA to the box center
        set_box (simu_struct.box, topol_struct, coord_vel_force_struct);

        if (crowder_struct.N_crwd [0] > 0) {
            printf ("Adding   %d  ions ...\n", crowder_struct.N_crwd [0]);
            //// Adding crowders
            int index = topol_struct.Nres_biomol;
            //        for (int i = 0; i < N_CRWD; i++) {
            for (int i = 0; i < 1; i++) {
                /// Calculate mean distance of crowders corresponding with the bulk concentration
                double mean_dist = cbrt (simu_struct.volume / crowder_struct.N_crwd [i]);
                // printf ("mean distance   %f\n", mean_dist);

                for (int j = 1; j <= crowder_struct.N_crwd [i]; j++) {
                    index++;
                    add_crowder (index, simu_struct.box, topol_struct, coord_vel_force_struct, 0, 0.8*mean_dist);
                }
            }
        }

        //////////////// Generate velocity ////////////////////
        for (int i = 1; i <= topol_struct.Nres; i++) {
            int Nsite;
            if (i > topol_struct.Nres_biomol) Nsite = topol_struct.Nsite [i];
            else Nsite = 2;
            for (int j = 0; j <= Nsite; j++)
                generate_atom_velocity (i, j, simu_struct.temp, topol_struct, coord_vel_force_struct);
        }

        Write_pdb_cg2 (mol_cg, topol_struct, crowder_struct, coord_vel_force_struct, simu_struct.ion_type, out_pdb);
    ///////////////////// RESTART //////////////////
    } else {
        read_rst_file (file_rst, topol_struct.Natm, simu_struct.box, coord_vel_force_struct);
        read_progress (file_info, simu_struct.progress, simu_struct.step_traj);
    }

    gen_con_matrix (topol_struct, mol_cg);
    gen_old_coord (topol_struct.Natm, simu_struct.dt, coord_vel_force_struct);

    simu_struct.R1_LIST = simu_struct.D_CUTOFF + simu_struct.dr1;
    simu_struct.R2_LIST = simu_struct.R1_LIST * simu_struct.R1_LIST;

/*    simu_struct.R1_LIST [0] = simu_struct.D_SR_CUTOFF + simu_struct.dr1;
    simu_struct.R2_LIST [0] = simu_struct.R1_LIST [0] * simu_struct.R1_LIST [0];

    simu_struct.R1_LIST [1] = simu_struct.D_CT_CUTOFF + simu_struct.dr1;
    simu_struct.R2_LIST [1] = simu_struct.R1_LIST [1] * simu_struct.R1_LIST [1];*/

//    tabulate_CT (EW_struct);
//    prepare_fourier_space (simu_struct, topol_struct.Natm, EW_struct);

    FILE * f_out;//, * f_dcd;
//    f_dcd = fopen (file_dcd, "w");
//    write_dcd_header (f_dcd, simu_struct, topol_struct.Natm);

//    dcdhandle *dcd;
    void *dcd = open_dcd_write (file_dcd, "dcd", simu_struct, topol_struct.Natm);
//    dcd = (dcdhandle *) v;

    f_out = fopen (file_output, "w");
//    FILE * f_error;
//    f_error = fopen (file_error, "w");

    double test = simu_struct.D_CUTOFF - 0.01;
    printf ("Electrostatic force  at %6.1f A: %16.8e\n", simu_struct.D_CUTOFF, DEBYE_HUCKEL_FORCE (test, simu_struct.kappaD, topol_struct.E_DH [1], simu_struct.D_CUTOFF, energy_struct.E_Q));
    printf ("Electrostatic energy at %6.1f A: %16.8e\n", simu_struct.D_CUTOFF, energy_struct.E_Q);

    if (simu_struct.progress == 0 && simu_struct.fix_solute == 0) {
        populate_lists (simu_struct, topol_struct, coord_vel_force_struct);
        update_list_hydbond (simu_struct, topol_struct, hydbond_struct);

//        prepare_kosinus_sinus (topol_struct, EW_struct, coord_vel_force_struct);
        deterministic_forces (mol_cg.res, simu_struct, topol_struct, hydbond_struct, stack_struct, coord_vel_force_struct, energy_struct, &myseed);

        fprintf (f_out, "Step      0\n");
//        fprintf (f_out, "        E_BOND           %16.8e\n", energy_struct.E_BOND);
//        fprintf (f_out, "        E_ANGLE          %16.8e\n", energy_struct.E_ANGLE);
//        fprintf (f_out, "        E_EXCLUDED       %16.8e\n", energy_struct.E_EXCLUDED);
        fprintf (f_out, "        E_NON-NAT_HB     %16.8e\n", energy_struct.E2_HB);
        fprintf (f_out, "        E_NAT_HB         %16.8e\n", energy_struct.E3_HB);
        fprintf (f_out, "        E_NON-NAT_STACK  %16.8e\n", energy_struct.E2_STACK);
        fprintf (f_out, "        E_NAT_STACK      %16.8e\n", energy_struct.E3_STACK);
        fprintf (f_out, "        E_electrostatic  %16.8e\n", energy_struct.E_Q);
        energy_struct.E_potential = energy_struct.E_BOND + energy_struct.E_ANGLE + energy_struct.E_EXCLUDED + energy_struct.E2_HB + energy_struct.E3_HB + energy_struct.E2_STACK + energy_struct.E3_STACK + energy_struct.E_Q;
        fprintf (f_out, "        E_potential      %16.8e\n", energy_struct.E_potential);

        double Rg = radius_gyration (topol_struct, coord_vel_force_struct);
        fprintf (f_out, "        Rg    %16.8e\n", Rg);

        if (topol_struct.tertiary_stack.size() > 0) {
            fprintf (f_out, "   Tertiary stacks considered\n");
            for (size_t i = 0; i < topol_struct.tertiary_stack.size(); i++) {
                fprintf (f_out, "%6d:::%4d    ", topol_struct.tertiary_stack[i].first, topol_struct.tertiary_stack[i].second);
                if (i % 10 == 9) fprintf (f_out, "\n");
            }
            fprintf (f_out, "\n");
        }

        fprintf (f_out, "\n");
        fflush (f_out);
//        write_dcd (coord_vel_force_struct, f_dcd, topol_struct.Natm);

        if (simu_struct.debug)
            debug_force (mol_cg.res, simu_struct, topol_struct, hydbond_struct, stack_struct, coord_vel_force_struct, energy_struct);

        // Coordinated-restraint MD
        if (simu_struct.restraint_step > 0) {
            printf ("Performing a quick fixed MD ...\n");
            for (int step = 1; step <= simu_struct.restraint_step; step++) {
//                prepare_kosinus_sinus (topol_struct, EW_struct, coord_vel_force_struct);
//                forces_PRMD (simu_struct, topol_struct, coord_vel_force_struct, EW_struct);

                for (int i = 1; i <= topol_struct.Natm; i++) {
                    coord_vel_force_struct.forcex [i] = 0;
                    coord_vel_force_struct.forcey [i] = 0;
                    coord_vel_force_struct.forcez [i] = 0;
                }
                non_bonded_interaction (simu_struct, topol_struct, coord_vel_force_struct, energy_struct);
                full_forces (topol_struct, coord_vel_force_struct, &myseed);
                move_crowders (simu_struct, topol_struct, coord_vel_force_struct);
//                simu_struct.relist_step++;
                check_shifts (simu_struct, topol_struct, hydbond_struct, coord_vel_force_struct);
            }
        }

    } else if (simu_struct.fix_solute && simu_struct.debug) {
        populate_lists (simu_struct, topol_struct, coord_vel_force_struct);
//        printf ("total interaction  %ld\n", simu_struct.neighborlist_mass);
//        prepare_kosinus_sinus (topol_struct, EW_struct, coord_vel_force_struct);
//        forces_PRMD (simu_struct, topol_struct, coord_vel_force_struct, EW_struct);

/*        for (int i = 1; i <= topol_struct.Natm; i++) {
            coord_vel_force_struct.forcex [i] = 0;
            coord_vel_force_struct.forcey [i] = 0;
            coord_vel_force_struct.forcez [i] = 0;
        }*/
        non_bonded_interaction (simu_struct, topol_struct, coord_vel_force_struct, energy_struct);
        debug_force_fixed (topol_struct.Natm, coord_vel_force_struct);
    }

    _Dcoordinate ext_force;
    if (simu_struct.pull_force) {
        double norm = sqrt (simu_struct.f_direction.x*simu_struct.f_direction.x + simu_struct.f_direction.y*simu_struct.f_direction.y + simu_struct.f_direction.z*simu_struct.f_direction.z);
        ext_force.x = simu_struct.f_direction.x / norm * simu_struct.force * PN_TO_KCALMOL;
        ext_force.y = simu_struct.f_direction.y / norm * simu_struct.force * PN_TO_KCALMOL;
        ext_force.z = simu_struct.f_direction.z / norm * simu_struct.force * PN_TO_KCALMOL;
    }

    populate_lists (simu_struct, topol_struct, coord_vel_force_struct);
    if (simu_struct.fix_solute == 0) update_list_hydbond (simu_struct, topol_struct, hydbond_struct);
//    simu_struct.relist_step = 0;

    double wall_time0 = get_wall_time ();
    double wall_time1 = wall_time0;

    printf ("Simulation started ...\n");

    //////// Regular simulation /////////
    if (simu_struct.fix_solute == 0)
        for (long int step = 1; step <= simu_struct.nstep; step++) {
//            prepare_kosinus_sinus (topol_struct, EW_struct, coord_vel_force_struct);
            deterministic_forces (mol_cg.res, simu_struct, topol_struct, hydbond_struct, stack_struct, coord_vel_force_struct, energy_struct, &myseed);
            full_forces (topol_struct, coord_vel_force_struct, &myseed);
            if (simu_struct.pull_force)
                apply_external_force (simu_struct.f_id1, simu_struct.f_id2, ext_force, coord_vel_force_struct);

            move_rigid_units (step, simu_struct, topol_struct, coord_vel_force_struct);
//            simu_struct.relist_step ++;
            simu_struct.progress ++;
            check_shifts (simu_struct, topol_struct, hydbond_struct, coord_vel_force_struct);

            if (step % simu_struct.step_traj == 0) {
                write_timestep (dcd, simu_struct.box, coord_vel_force_struct);
                write_rst (file_rst, topol_struct.Natm, simu_struct.box, coord_vel_force_struct);
            }

            if (step % simu_struct.step_energy == 0) {
                double Rg = radius_gyration (topol_struct, coord_vel_force_struct);
                fprintf (f_out, "Step    %lld\n", simu_struct.progress);
//                fprintf (f_out, "        E_BOND           %16.8e\n", energy_struct.E_BOND);
//                fprintf (f_out, "        E_ANGLE          %16.8e\n", energy_struct.E_ANGLE);
//                fprintf (f_out, "        E_EXCLUDED       %16.8e\n", energy_struct.E_EXCLUDED);
                fprintf (f_out, "        E_NON-NAT_HB     %16.8e\n", energy_struct.E2_HB);
                fprintf (f_out, "        E_NAT_HB         %16.8e\n", energy_struct.E3_HB);
                fprintf (f_out, "        E_NON-NAT_STACK  %16.8e\n", energy_struct.E2_STACK);
                fprintf (f_out, "        E_NAT_STACK      %16.8e\n", energy_struct.E3_STACK);
                fprintf (f_out, "        E_electrostatic  %16.8e\n", energy_struct.E_Q);
                energy_struct.E_potential = energy_struct.E_BOND + energy_struct.E_ANGLE + energy_struct.E_EXCLUDED + energy_struct.E2_HB + energy_struct.E3_HB + energy_struct.E2_STACK + energy_struct.E3_STACK + energy_struct.E_Q;
                fprintf (f_out, "        E_potential      %16.8e\n", energy_struct.E_potential);

                fprintf (f_out, "        Rg     %16.8e  A\n", Rg);

                if (stack_struct.s3_N > 0) {
                    fprintf (f_out, "    Individual tertiary stack energy:\n");

                    for (int i = 1; i <= stack_struct.s3_N; i++) {
                        fprintf (f_out, "%+7.3f   ", stack_struct.st_energy [i]);
                        if (i % 10 == 0) fprintf (f_out, "\n");
                    }
                    fprintf (f_out, "\n");
                }

                fprintf (f_out, "\n");
                fflush (f_out);

//                printf ("Temperature   %6.1f\n", calc_temp (topol_struct, coord_vel_force_struct));
                double wall_time2 = get_wall_time ();
                double dtime1 = wall_time2 - wall_time0;
                double dtime2 = wall_time2 - wall_time1;
                write_info (file_info, step, dtime1, dtime2, topol_struct.bead_zeroQ, Rg, simu_struct, crowder_struct, energy_struct);
                wall_time1 = wall_time2;
            }
        }
    else
        /////// Fixed-solute simulation ///////////
        for (long int step = 1; step <= simu_struct.nstep; step++) {
//            prepare_kosinus_sinus (topol_struct, EW_struct, coord_vel_force_struct);
//            forces_PRMD (simu_struct, topol_struct, coord_vel_force_struct, EW_struct);

            for (int i = 1; i <= topol_struct.Natm; i++) {
                coord_vel_force_struct.forcex [i] = 0;
                coord_vel_force_struct.forcey [i] = 0;
                coord_vel_force_struct.forcez [i] = 0;
            }
            non_bonded_interaction (simu_struct, topol_struct, coord_vel_force_struct, energy_struct);
            full_forces (topol_struct, coord_vel_force_struct, &myseed);
            move_crowders (simu_struct, topol_struct, coord_vel_force_struct);
//            simu_struct.relist_step++;
            simu_struct.progress++;
            check_shifts (simu_struct, topol_struct, hydbond_struct, coord_vel_force_struct);

            if (step % simu_struct.step_traj == 0) {
                write_timestep (dcd, simu_struct.box, coord_vel_force_struct);
                write_rst (file_rst, topol_struct.Natm, simu_struct.box, coord_vel_force_struct);
            }

            if (step % simu_struct.step_energy == 0) {
                double wall_time2 = get_wall_time ();
                double dtime1 = wall_time2 - wall_time0;
                double dtime2 = wall_time2 - wall_time1;
                double Rg = 0.;
                write_info (file_info, step, dtime1, dtime2, topol_struct.bead_zeroQ, Rg, simu_struct, crowder_struct, energy_struct);
                wall_time1 = wall_time2;
            }
        }

//    fclose (f_error);
    fclose (f_out);
//    fclose (f_dcd);
    close_file_write (dcd);

//    for (int i = 1; i <= 10; i++) free (simu_struct.pair_potential [i]);
//    free (simu_struct.pair_potential);
    for (int i = 1; i <= topol_struct.Nres_biomol; i++) free (topol_struct.atom_key [i]);
    free (topol_struct.atom_key);

    free_stacks (stack_struct);
    free_hydrogen_bonds (topol_struct, hydbond_struct);

   // free (coord_vel_force_struct.coordx);
   // free (coord_vel_force_struct.coordy);
   // free (coord_vel_force_struct.coordz);
   // free (coord_vel_force_struct.coord_oldx);
   // free (coord_vel_force_struct.coord_oldy);
   // free (coord_vel_force_struct.coord_oldz);
   // free (coord_vel_force_struct.velox);
   // free (coord_vel_force_struct.veloy);
   // free (coord_vel_force_struct.veloz);
   // free (coord_vel_force_struct.forcex);
   // free (coord_vel_force_struct.forcey);
   // free (coord_vel_force_struct.forcez);
//    free (coord_vel_force_struct.flx);
//    free (coord_vel_force_struct.fly);
//    free (coord_vel_force_struct.flz);

/*    for (int i = 0; i <= EW_struct.N_vec; i++) {
        free (EW_struct.sinus   [i]);
        free (EW_struct.kosinus [i]);
    }
    free (EW_struct.sinus);         free (EW_struct.kosinus);*/

    for (int i = 1; i <= simu_struct.Mc_total; i++)
        free (simu_struct.cell_content [i]);
    free (simu_struct.cell_content); //free (simu_struct.cell_mass);
    return (0);
}
