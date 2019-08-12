#ifndef ANAL_H
#define ANAL_H

double radius_gyration (_topol_struct                 &topol_struct,
                        const _coord_vel_force_struct &coord_vel_force_struct) {
    _Dcoordinate COM;
    COM.x = 0;       COM.y = 0;        COM.z = 0;

    int Nbead = topol_struct.Natm_biomol;

    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        if (i == topol_struct.ligand_bead) {
            Nbead--;
            continue;
        }

        COM.x += coord_vel_force_struct.coordx[i];
        COM.y += coord_vel_force_struct.coordy[i];
        COM.z += coord_vel_force_struct.coordz[i];
    }
    COM.x /= Nbead;     COM.y /= Nbead;     COM.z /= Nbead;
    double x = coord_vel_force_struct.coordx [topol_struct.bead_center] - COM.x;
    double y = coord_vel_force_struct.coordy [topol_struct.bead_center] - COM.y;
    double z = coord_vel_force_struct.coordz [topol_struct.bead_center] - COM.z;
    double d2 = x*x + y*y + z*z;

    double Rg = 0;
    double min = 9e9;
    int tmp = 0;
    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        if (i == topol_struct.ligand_bead) continue;

        double dx = COM.x - coord_vel_force_struct.coordx[i];
        double dy = COM.y - coord_vel_force_struct.coordy[i];
        double dz = COM.z - coord_vel_force_struct.coordz[i];
        double d = dx*dx + dy*dy + dz*dz;
        Rg += d;
        if (d < min && sqrt(d2) - sqrt(d) >= 10.) {
            min = d;
            tmp = i;
        }
    }
    if (tmp != 0) topol_struct.bead_center = tmp;
    return sqrt (Rg/Nbead);
}

#endif
