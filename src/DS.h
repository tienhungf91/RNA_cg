#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

#include <string>
#include <vector>
#include "const.h"

struct _Dcoordinate {
//    double x = 0, y = 0, z = 0;
    double x, y, z;
};

struct _Icoordinate {
//    int x = 0, y = 0, z = 0;
    int x, y, z;
};

struct _Atom {
    size_t index;
    std::string name;
    std::string element;
    double mass;
    _Dcoordinate coord;
};

struct _Res_aa {
    std::string name;
    std::vector<_Atom> atom;
    size_t index;
    std::string chain;
    std::vector<size_t> connectedTo;
};

struct _Mol_aa {
    std::vector<_Res_aa> res;
    std::string type;
};

struct _Bead {
    std::string name;
//    std::string resname;
    size_t index;
//    size_t res_index;
    std::string chain;
    double mass;
    double charge;
	std::vector<size_t> connectedTo;
    _Dcoordinate coord;
};

struct _Res_cg {
    std::string name;
    size_t index;
    std::vector<_Bead> bead;
    std::string chain;
    double mass;
    double charge;
    std::vector<size_t> connectedTo;
};

struct _Mol_cg {
    std::vector<_Res_cg> res;
    std::string type;
    int netcharge;
};

struct _dx {
    std::string type;
    double conc;
    std::vector<double> value;
    _Icoordinate ngrid;
    _Dcoordinate origin;
    _Dcoordinate delta;
};

struct _simu_struct {
    int step_traj;
    int step_energy;
    long int nstep;
    long int restraint_step;
    int step_adj_box;
    
    //Box parameters
    double box [3];
    double volume;

    //General simulation parameters
    double visc_gamma;
    double dt;
//    double dtq2 = 0.5 * dt * dt;
    double temp;
    long long progress;

    bool pull_force;
    double force;
    _Dcoordinate f_direction;
    int f_id1, f_id2;

    bool restart;  // Restart run
    bool debug;    // Debug run
    bool fix_solute;
    int seed;

    ///////////Electrostatics/////////
    bool Coulomb;
    double lB_T; //Bjerrum length lB = q^2 / kT
    double Ze; //scaled phosphate charge
    double epsilon, d_epsilon; // dielectric constant and d_ep/dT
    double kappaD, lambdaD;  // lambda = Debye length
                             // kappa  = 1./Debye length
    double b;   // Length per unit (bare) charge in the RNA (for charge scaling purpose)
    int ion_q, ion_type;

//    double D_SR_CUTOFF;
//    double D_CT_CUTOFF;
//    double D2_CT_CUTOFF;
    double D_CUTOFF, D2_CUTOFF;
    double lB_D2; // lB / D2_CUTOFF
//    int LR_skip;

    ///////////// Cells (for list population or pcfs calculation)//////////////
    _Icoordinate Mc;
    int Mc_total; //cell numbers

    std::vector<int> cell_mass;
//    std::vector< std::vector<int> > cell_content;
    int ** cell_content;//, * cell_mass = NULL; //cell content and mass

    /////////////////////////////////Lists////////////////////////////////////
//    int * list_content [NPI];
//    int list_mass [NPI];
//    int list_mass_SR, list_mass_LR;
    long int HBneighborlist_mass, DHneighborlist_mass;
//    int * list_content_vdw, * list_content_ele;
//    std::vector<int> list_content_SR, list_content_LR;
    std::vector<long int> HBneighborlist, DHneighborlist;

//    double R1_LIST [NPI];
//    double R2_LIST [NPI];
    double R1_LIST, R2_LIST;

    double dr1;
    double relist;
//    int relist_step;
    std::vector<double> pair_potential;
};

///////////////////////////////////////////
struct _coord_vel_force_struct {
    //Coordinates, velocities and forces
/*    std::vector<_Dcoordinate> coord;
    std::vector<_Dcoordinate> coord_old;
    std::vector<_Dcoordinate> velo;
    std::vector<_Dcoordinate> force;
    std::vector<_Dcoordinate> fl;*/
//    double * coordx;
//    double * coordy;
//    double * coordz;
    std::vector<double> coordx, coordy, coordz;

//    double * coord_oldx;
//    double * coord_oldy;
//    double * coord_oldz;
    std::vector<double> coord_oldx, coord_oldy, coord_oldz;

//    double * velox;
//    double * veloy;
//    double * veloz;
    std::vector<double> velox, veloy, veloz;

    double * forcex;
    double * forcey;
    double * forcez;

    double * flx;
    double * fly;
    double * flz;

    std::vector<double> maxwell_force;
};

////////////////////////////////////////////////
//////// Index starts at 1 for RNA
struct _topol_struct {
    //Particle numbers
    int N_amino;       // number of residues in the protein
    int Nres_biomol;    // number of backbone atoms in the protein 
    int Nres;          // total number of backbone atoms in the protein and crowders
    std::vector<int> Nsite;
//    int * Nsite = NULL;   // number of sites per backbone atom or crowder ("Nsite = m_crwd - 1" for crowders)
    int Natm_biomol;   // total number of atoms in RNA
    int Npair_biomol;   // number of different pairs of atoms in RNA
    int Natm;         // total number of atoms in RNA and crowders
    int ligand_index, ligand_bead;

    int bead_center;
//    std::vector< std::vector<int> > atom_key;
    std::vector<int> bead_zeroQ;
    std::vector<int> res_pdb2res_simu;             // map pdb res -> simulation res
    std::vector<int> res_simu2res_pdb;             // map simulation res -> pdb res
    std::vector<int> stack_ignore;
    std::vector<int> hbond_ignore;
    std::vector<std::pair<int, int> > secondary_stack, tertiary_stack;
    int ** atom_key, ** con_matrix;
//    std::vector<char> amino_key;
    std::vector<int> ter_card;         // Specify where the TER card is
//    char * amino_key = NULL;
    std::vector<int> part_key, maxi_key;// INDX, JNDX;
//    int * part_key = NULL;
//    int * maxi_key = NULL;
//    int * INDX = NULL, * JNDX = NULL;
//    std::vector<int> residue_mass;
//    int * residue_mass = NULL;  /// 1D array contains how many A, G, C, U are
                                 // 0 1 2 3 4 -> P A G C U
//    std::vector< std::vector<int> > residue_content;
//    int ** residue_content = NULL;  /// 2D array contains indexs of residue type
                                 // 0 1 2 3 4 -> P A G C U

//    std::vector<int> RIGID_SET, RIGID_END;
//    int * RIGID_SET = NULL, * RIGID_END = NULL;

    //Individual properties of atoms
//    double D [NS] = { 4.2, 5.8, 6.0, 1.5852, 3.8960, 5.3160 }; //phosphate, RIBOSE, BASE, Mg, Cl, K
//    double PD [NPS];
//    double PD_sq [NPS];

    /////////////////////////////////////
    double MASS [NA + 1]; //atom masses
    double RLJ [NA + 1]; //atom LJ radii
    double ELJ [NA + 1]; //atom LJ well depths
    double QC [NA + 1]; //atom charges
    double VISC [NA + 1];
    double K1 [NA + 1];
    double K2 [NA + 1];
    double SIGMA_FORCE [NA + 1];
    double D_LJ [NPA + 1]; //LJ pair diameters
    double D2_LJ [NPA + 1]; //LJ square pair diameters
    double D2_LJ_CUTOFF [NPA + 1]; //LJ square cutoff distances
    double D2_LJ_OVERLAP [NPA + 1]; //LJ square overlap distances
    double D2_LJ_MINIMUM [NPA + 1]; //LJ square minimum distances
    double E_LJ [NPA + 1]; //LJ energy prefactors
    double F_LJ [NPA + 1]; //LJ force prefactors
    double Q_C [NPA + 1]; //Coulomb prefactors
    double E_DH [NPA + 1];

/*    std::vector<double> CMS, VCMS, RMASS, RVISC, RXYZ;
    std::vector<double> AXES, W, IR1, IR2, IR3;
    std::vector<double> AV, BV, CV, FV, GV, HV;
    std::vector< std::vector<double> > RX, RY, RZ;*/
/*    //Centre of mass motion
    double * CMS = NULL;
    double * VCMS = NULL;
    double * RMASS = NULL;
    double * RVISC = NULL;
    double * RXYZ = NULL;

    //Rotational motion
    double * AXES = NULL;
    double * W = NULL;
    double ** RX = NULL, ** RY = NULL, ** RZ = NULL;
    double * IR1 = NULL, * IR2 = NULL, * IR3 = NULL;
    double * AV = NULL, * BV = NULL, * CV = NULL;
    double * FV = NULL, * GV = NULL, * HV = NULL;*/
};

//////////////////////////////////////////////////////
///////// Index starts at 0 for crowders
struct _crowder_struct {
    //Crowders
    double conc [N_CRWD]; //molar concentrations of crowders
//    double phi_crwd [N_CRWD] = { 0.0, 0.0, 0.0 }; //volume fractions of crowders
//    double vol [N_CRWD]; //crowder volumes
    std::vector<int> crowder_key;
//    int * crowder_key = NULL;
    int N_crwd [N_CRWD]; //the number of crowders of each type
    int Nsite [N_CRWD]; //the number of sites in a crowder of each type

    std::vector<int> crowder_mass;
    std::vector< std::vector<int> > crowder_content;
//    int * crowder_mass = NULL;
//    int ** crowder_content = NULL;

//    double * RX_crwd [N_CRWD]; //x coordinates of sites in a crowder of each type
//    double * RY_crwd [N_CRWD]; //y coordinates of sites in a crowder of each type
//    double * RZ_crwd [N_CRWD]; //z coordinates of sites in a crowder of each type
    /// It seems that the code does not use RX_crwd, RY_crwd and RZ_crwd to store crowder coords,
    /// instead, they are stored in x, y and z in coord_vel_force_struct

//    int * part_key_crwd [N_CRWD]; //site (atom) "part_key" in a crowder of each type
//    int * maxi_key_crwd [N_CRWD]; //site (atom) "maxi_key" in a crowder of each type
    std::vector< std::vector<int> > part_key_crwd, maxi_key_crwd;
//    int RIGID_SET_crwd [N_CRWD]; //the first site in a crowder rigid sub-unit
//    int RIGID_END_crwd [N_CRWD]; //the last site in a crowder rigid sub-unit
};

////////////////////////////////////////////////////
struct _stack_struct {
    bool eval;
    //Stacks
    double st_D; //adjustment for tertiary stacks
    double ss_D; //adjustment for secondary stacks
    int s3_N; //number of tertiary stacks
    int st_N; //total number of stacks

    double score_cutoff;
/*    int * st_i = NULL;
    int * st_j = NULL;

    double * st_r = NULL;
    double * st_theta1 = NULL;
    double * st_theta2 = NULL;
    double * st_psi = NULL;
    double * st_psi1 = NULL;
    double * st_psi2 = NULL;
    double * st_E = NULL;
    double * st_psi0 = NULL;
    double * st_psi10 = NULL;
    double * st_psi20 = NULL;
    double * st_r2 = NULL;

    double * st_energy = NULL;
    int * st_status = NULL;
    int * ST_EXCESS = NULL;
    int * ST_ATOM_N = NULL;*/

    std::vector<double> st_r, st_theta1, st_theta2, st_psi, st_psi1, st_psi2, st_psi0, st_psi10, st_psi20;
    std::vector<double> st_E, st_energy, e3_stack;
    std::vector<int> st_i3, st_j3, st_i2, st_j2, st_status, ST_ATOM_N, Nstack_allowed, Nstack_formed;
    std::vector<bool> ligand_stack;
//    int ** ST_ATOM;
//    std::vector< std::vector<int> > ST_ATOM;
//    int ** ATOM_ST = NULL;
//    std::vector< std::pair<int, int> > ATOM_ST;
//    double * e3_stack = NULL;
    double k_r2, k_r3, k_theta3, k_psi2, k_psi3;
};

////////////////////////////////////////////////////////
struct _hydbond_struct {
    bool eval;
    //Hydrogen bonds
    double hs_D;
    int hb_N;
    std::vector<int> hb_dode, hb_k10, hb_k11, hb_k12, hb_k20, hb_k21, hb_k22;
/*    int * hb_dode = NULL;
    int * hb_k10 = NULL;
    int * hb_k11 = NULL;
    int * hb_k12 = NULL;
    int * hb_k20 = NULL;
    int * hb_k21 = NULL;
    int * hb_k22 = NULL;

    double * hb_r = NULL;
    double * hb_theta1 = NULL;
    double * hb_theta2 = NULL;
    double * hb_psi = NULL;
    double * hb_psi1 = NULL;
    double * hb_psi2 = NULL;
    double * hb_E = NULL;
    double * hb_psi0 = NULL;
    double * hb_psi10 = NULL;
    double * hb_psi20= NULL;
    double * hb_energy = NULL;
    int * hb_status = NULL;*/

    std::vector<double> hb_r, hb_theta1, hb_theta2, hb_psi, hb_psi1, hb_psi2;
    std::vector<double> hb_E, hb_psi0, hb_psi10, hb_psi20, hb_energy, hb_f;
    std::vector<int> hb_status;
    std::vector< std::pair <int, int> > resid;
    std::vector<bool> ligand_bond;

//    std::vector<int> hb_RES1, hb_RES2, HB_Nsite; //HB_INDX, HB_JNDX,
    std::vector<int> HB_Nsite;
    std::vector<int> VALENCE, EXCESS, HB_EXCESS;
//    int * hb_RES1 = NULL;
//    int * hb_RES2 = NULL;
    char ** hb_code;
//    int * HB_Nsite = NULL; //the number of hydrogen bonding atoms in a "sugar + base" or a phosphate
    int HB_Natm;      //the total number of hydrogen bonding atoms
//    std::vector< std::vector<int> > HB_atom_key;
    int ** HB_atom_key;
//    int * HB_INDX = NULL, * HB_JNDX = NULL;
//    int * VALENCE = NULL;
//    int * HB_EXCESS = NULL;
//    int * EXCESS = NULL;
    double weight_nonnative;

//    int * HB_ATOM_N = NULL;
    std::vector<int> HB_ATOM_N, ATOM_HB_N, HB_PAIR_N, HB_A, hb_K;
    int ** HB_ATOM;
//    std::vector< std::vector<int> > HB_ATOM, ATOM_HB;
//    int * ATOM_HB_N = NULL;
    int ** ATOM_HB;
//    int * HB_PAIR_N = NULL;
    int ** HB_PAIR;
//    std::vector< std::vector<int> > HB_PAIR;
//    int * HB_A = NULL;
//    int * hb_K = NULL;
    double k_r, k_theta, k_phi;
};

///////////////////////////////////////////////////////////
struct _energy_struct {
    //Energies
    double E_BOND;
    double E_ANGLE;
    double E_EXCLUDED;
    double E2_HB;    // non-native bond energy
    double E3_HB;    // native     bond energy

    double E2_STACK; // secondary stacks
    double E3_STACK; // tertiary stacks
    double E_Q;      // Electrostatic
    double E_potential;
/*    double ELJ_NON_NATIVE = 0;
    double E_YUKAWA = 0;
    double E_BOND_B = 0;
    double E_BOND_CRWD = 0;
    double E_VALENCE_B = 0;
    double E_VALENCE_CRWD = 0;
    double E_DIHEDRAL_B = 0;
    double E_DIHEDRAL_CRWD = 0;*/
};

/////////////////////////////////////////////////////////////
struct _EW_struct {
    /////////Real Space Ewald/////////
    double sigma, kappa; //sigma < side_X / 5
    int CTM;
    std::vector<double> CTE_array, CTF_array, sin_array, cos_array;

    //////Fourier Space Ewald/////////
    int N_vec;
    std::vector<double> CONST_EWALD, fourier_x, fourier_y, fourier_z;
    std::vector<double> A1, A2;
    double ** sinus, ** kosinus;
    double EXP_CUTOFF, CONST_ADD, CONST_REMOVE;
    //Chemical potential
//    double CHEMICAL_POTENTIAL = 0;
//    long int COUNT_CH = 0;
};

#endif
