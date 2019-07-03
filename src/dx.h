#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include "DS.h"

// ********************************************
// ***            Read dx file              ***
// ********************************************
void read_dx (std::string &dx_file,
              _dx &dx) {

    double tmp1, tmp2, tmp3;
    std::ifstream DXFILE (dx_file.c_str());
    if (DXFILE.is_open()) {
        std::string line;
        while (getline(DXFILE, line)) {
            std::istringstream iss(line);
            while (iss >> tmp1)
                dx.value.push_back(tmp1);
            if (line.find("object 1 class gridpositions counts") != std::string::npos) {
                line.replace(0,35," ");
                std::istringstream iss(line);
//                dx.ngrid.resize(3);
                iss >> dx.ngrid.x >> dx.ngrid.y >> dx.ngrid.z;
            } else if (line.find("origin") != std::string::npos) {
                line.replace(0,6," ");
                std::istringstream iss(line);
//                dx.origin.resize(3);
                iss >> dx.origin.x >> dx.origin.y >> dx.origin.z;
            } else if (line.find("delta") != std::string::npos) {
                line.replace(0,5," ");
                std::istringstream iss(line);
//                dx.delta.resize(3);
                iss >> tmp1 >> tmp2 >> tmp3;
                if (tmp1 != 0)
                    dx.delta.x = tmp1;
                else if (tmp2 != 0)
                    dx.delta.y = tmp2;
                else if (tmp3 != 0)
                    dx.delta.z = tmp3;
            }
        }
        DXFILE.close();
        if (dx.ngrid.x * dx.ngrid.y * dx.ngrid.z != dx.value.size()) {
            std::cout << "Number of grid points in " << dx_file << " not matched!!\n";
            std::cout << "Number of grid points (x y z):\t" << dx.ngrid.x << "\t" << dx.ngrid.y << "\t" << dx.ngrid.z << std::endl;
            std::cout << "Number of values:\t" << dx.value.size() << std::endl;
            exit (0);
        }
    } else {
        std::cerr << "Unable to open file " << dx_file << std::endl;
        exit (0);
    }
}

// ********************************************
// ***           Write dx file              ***
// ********************************************

void write_dx (const _dx &dx,
               std::string &outfile) {

    std::ofstream OUT (outfile.c_str());
    if (OUT.is_open()) {
        OUT << "object 1 class gridpositions counts" << std::setw(8) << dx.ngrid.x << std::setw(8) << dx.ngrid.y << std::setw(8) << dx.ngrid.z << std::endl;
        OUT << "origin " << std::setw(15) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.origin.x <<\
                            std::setw(15) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.origin.y <<\
                            std::setw(15) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.origin.z << std::endl;

        OUT << "delta " << std::setw(16) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.delta.x << " 0 0\n";
        OUT << "delta  0" << std::setw(16) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.delta.y << " 0\n";
        OUT << "delta  0 0" << std::setw(16) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.delta.z << std::endl;
        OUT << "object 2 class gridconnections counts" << std::setw(8) << dx.ngrid.x << std::setw(8) << dx.ngrid.y <<\
            std::setw(8) << dx.ngrid.z << std::endl;
        OUT << "object 3 class array type double rank 0 items " << dx.ngrid.x*dx.ngrid.y*dx.ngrid.z << " data follows";

        for (size_t i = 0; i < dx.value.size(); i++) {
            if (i % 3 == 0) OUT << std::endl;
            OUT << std::setw(16) << std::setprecision(5) << std::setiosflags(std::ios::fixed) << std::scientific << dx.value[i];
        }

        OUT << "\nobject \"Untitled\" call field";
        OUT.close();
    } else {
        std::cout << "Can't write to file " << outfile << std::endl;
        exit (0);
    }
}

/////////////////////////////////////////////
// Check to see dx grids are different or not
/////////////////////////////////////////////
/*bool check_dx (const std::vector<_dx> &dx) {
    bool check = 0;
    for (size_t i = 1; i < dx.size(); i++) {
        if (dx[i].ngrid != dx[0].ngrid) {
            check = 1;
            break;
        }
        if (dx[i].origin != dx[0].origin) {
            check = 1;
            break;
        }
        if (dx[i].delta != dx[0].delta) {
            check = 1;
            break;
        }
    }
    return check;
}*/

// *****************************************************************
// *****    Get grid point coordinate from the 1d index        *****
// *****************************************************************
_Dcoordinate dx_1Dindex2coord (const _dx &dx,
                               size_t index) {

    size_t x_index = floor (index/(dx.ngrid.y*dx.ngrid.z));
    size_t y_index = floor ((index - x_index * dx.ngrid.y * dx.ngrid.z) / dx.ngrid.z);
    size_t z_index = index - x_index*dx.ngrid.y*dx.ngrid.z - y_index*dx.ngrid.z;

    _Dcoordinate coord;
    coord.x = dx.origin.x + dx.delta.x * x_index;
    coord.y = dx.origin.y + dx.delta.y * y_index;
    coord.z = dx.origin.z + dx.delta.z * z_index;

    return coord;
}

// ***************************************************
// *****    Return the 1D index from 3D index  *******
// ***************************************************
long dx_3Dindex_to_1Dindex (const _dx &dx,
                            int x, int y, int z) {
    return x*dx.ngrid.y*dx.ngrid.z + y*dx.ngrid.z + z;
}

/////////////////////////////////////////////////////////
//         Return the 1D index from 3D coordinate      //
/////////////////////////////////////////////////////////
long dx_3Dcoord_to_1Dindex (const _dx &dx,
                            float &x, float &y, float &z) {
    int i = (int) ((x - dx.origin.x) / dx.delta.x);
    int j = (int) ((y - dx.origin.y) / dx.delta.y);
    int k = (int) ((z - dx.origin.z) / dx.delta.z);
    return dx_3Dindex_to_1Dindex (dx, i, j, k);
}

// ***************************************************************
// ***      Extract a smaller dx from the original dx          ***
// ***************************************************************
void extract_dx (const _dx &dx, _dx &dx_out,
                 _Dcoordinate &min,
                 _Dcoordinate &max) {
    dx_out.delta = dx.delta;

    // Get the small and large indexes
    size_t imin = ceil ((min.x - dx.origin.x) / dx.delta.x);
    size_t jmin = ceil ((min.y - dx.origin.y) / dx.delta.y);
    size_t kmin = ceil ((min.z - dx.origin.z) / dx.delta.z);

//    dx_out.origin.resize (3);
    dx_out.origin.x = dx.origin.x + imin*dx.delta.x;
    dx_out.origin.y = dx.origin.y + jmin*dx.delta.y;
    dx_out.origin.z = dx.origin.z + kmin*dx.delta.z;

    size_t imax = floor ((max.x - dx.origin.x) / dx.delta.x);
    size_t jmax = floor ((max.y - dx.origin.y) / dx.delta.y);
    size_t kmax = floor ((max.z - dx.origin.z) / dx.delta.z);

//    dx_out.ngrid.resize (3);
    dx_out.ngrid.x = imax - imin + 1;
    dx_out.ngrid.y = jmax - jmin + 1;
    dx_out.ngrid.z = kmax - kmin + 1;

    // Loop between the small and large indexes to get value for the extracted dx
    for (size_t i = imin; i <= imax; i++)
        for (size_t j = jmin; j <= jmax; j++)
            for (size_t k = kmin; k <= kmax; k++) {
                size_t index = (size_t) (i*dx.ngrid.y*dx.ngrid.z + j*dx.ngrid.z + k);
                dx_out.value.push_back (dx.value[index]);
            }

    // Sanity check
    if (dx_out.ngrid.x * dx_out.ngrid.y * dx_out.ngrid.z != dx_out.value.size()) {
        std::cerr << "Number of grid point not matched\n";
        exit (0);
    }
}

// *****************************************************************
// ****         Adjust values of the dx: new = a*old + b         ***
// *****************************************************************
_dx scale_dx (const _dx &dx,
              double a,
              double b) {
    _dx result;
    result.value.resize (dx.value.size());
    result.ngrid = dx.ngrid;
    result.origin = dx.origin;
    result.delta = dx.delta;

    for (size_t i = 0; i < dx.value.size(); i++)
        result.value[i] = a*dx.value[i] + b;

    return result;
}

// *************************************************************
// ****        Integrate all grid points of a dx       *********
// *************************************************************
double integral_dx (const _dx &dx) {
    double result = 0;
    for (size_t i = 0; i < dx.value.size(); i++)
        result += dx.value[i];
    return result;
}
