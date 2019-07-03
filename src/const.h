#ifndef CONST_H
#define CONST_H
const double PI        = 3.14159265358979323846264338328;
const double E         = 2.71828182845904523536028747135;
const double TWOPI     = 6.2831853071795865;
const double HALFPI    = 1.5707963267948966;
const double ONE_THIRD = 0.33333333333333333333333333333;
const double ENERGY_BLOW = 50.;
const double AVOGADRO = 6.0221409e23;
const int NSTEP_RELIST = 40;  // seems to work well with dt=0.05
const double KB = 0.001987204118;   // Boltzmann constant in kcal/mol.K
const double KELVIN_TO_KT = 1.9872057946353295439632e-3;  // Convert Kelvin temp to kT in kcal/mol
const double PN_TO_KCALMOL = 0.0143932620936902486;  // Convert force in pN to comply with energy in kcal/mol (with distance in Angstrom)

const int NPI = 2; //the number of different interactions which require a list
               //LJ (WCA) -> InInt = 0; Coulomb -> InInt = 1
               //used in topol_struct
const int NS = 6; ///// phosphate, RIBOSE, BASE, Mg, Cl, K
//const int NPS = NS * (NS + 1) / 2; //pairs of atoms
                       /// used in topol_struct
const int NA = 9; //the number of different "amber" atom types, which is different from ns above
                  // because this distinguish different base types
                  // and which is the number of entries in MAXIKEY file/the number of different "amber" atom types
const int NPA = NA * (NA + 1) / 2; //pairs of "amber" atoms
                            ///// used in topol_struct
const int N_CRWD = 3; //the number of different crowder types
                      // used in crowder_struct
#endif
