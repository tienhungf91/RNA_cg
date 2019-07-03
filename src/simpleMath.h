#ifndef SMATH_H
#define SMATH_H

_Dcoordinate cross_product (const _Dcoordinate &v1,
                            const _Dcoordinate &v2) {
    _Dcoordinate res;
    res.x = v1.y * v2.z - v1.z * v2.y;
    res.y = v1.z * v2.x - v1.x * v2.z;
    res.z = v1.x * v2.y - v1.y * v2.x;
    return res;
}

/////////////////////////////////////////////////////
void cross_product2 (double * v1, double * v2, double * res) {
	res [0] = v1 [1] * v2 [2] - v1 [2] * v2 [1];
	res [1] = v1 [2] * v2 [0] - v1 [0] * v2 [2];
	res [2] = v1 [0] * v2 [1] - v1 [1] * v2 [0];
}

double norm_vector (const _Dcoordinate &v) {
    return sqrt (v.x*v.x + v.y*v.y + v.z*v.z);
}

////////////////////////////////////////////////////////
double dist_coord (const _Dcoordinate &a, const _Dcoordinate &b) {
    double x = a.x - b.x;
    double y = a.y - b.y;
    double z = a.z - b.z;
    return sqrt (x*x + y*y + z*z);
}

//////////////////////////////////////////////////////////
double compute_angle (const _Dcoordinate &v1, const _Dcoordinate &v2) {
    return 180./PI * acos ((v1.x*v2.x + v1.y*v2.y + v1.z*v2.z) / (norm_vector(v1)*norm_vector(v2)));
}

///////////////////////////////////////////////////////////////////////
_Dcoordinate gen_vector (const _Dcoordinate &a, const _Dcoordinate &b) {
    _Dcoordinate result;
    result.x = b.x - a.x;
    result.y = b.y - a.y;
    result.z = b.z - a.z;
    return result;
}

///////////////////////////////////////////////////////////
double CC_distance (int i1, int i2,
                    const _coord_vel_force_struct &coord_vel_force_struct) {
    double x = coord_vel_force_struct.coordx [i1] - coord_vel_force_struct.coordx [i2];
    double y = coord_vel_force_struct.coordy [i1] - coord_vel_force_struct.coordy [i2];
    double z = coord_vel_force_struct.coordz [i1] - coord_vel_force_struct.coordz [i2];
    return sqrt (x*x + y*y +z*z);
}

//////////////////////////////////////////////////
double valence_angle (int i1, int i2, int i3,
                      const _coord_vel_force_struct &coord_vel_force_struct) {
    double x1 = coord_vel_force_struct.coordx [i1] - coord_vel_force_struct.coordx [i2];
    double y1 = coord_vel_force_struct.coordy [i1] - coord_vel_force_struct.coordy [i2];
    double z1 = coord_vel_force_struct.coordz [i1] - coord_vel_force_struct.coordz [i2];

    double x2 = coord_vel_force_struct.coordx [i3] - coord_vel_force_struct.coordx [i2];
    double y2 = coord_vel_force_struct.coordy [i3] - coord_vel_force_struct.coordy [i2];
    double z2 = coord_vel_force_struct.coordz [i3] - coord_vel_force_struct.coordz [i2];

    double norm1 = x1*x1 + y1*y1 + z1*z1;
    double norm2 = x2*x2 + y2*y2 + z2*z2;

    double norm12 = sqrt (norm1 * norm2);

    double kosinus = (x1*x2 + y1*y2 + z1*z2) / norm12;
    return acos (kosinus);
}

/////////////////////////////////////////////////
double dihedral_angle (int i1, int i2, int i3, int i4,
                       const _coord_vel_force_struct &coord_vel_force_struct) {
    _Dcoordinate v1, v2, v3;
    v1.x = coord_vel_force_struct.coordx [i2] - coord_vel_force_struct.coordx [i1];
    v1.y = coord_vel_force_struct.coordy [i2] - coord_vel_force_struct.coordy [i1];
    v1.z = coord_vel_force_struct.coordz [i2] - coord_vel_force_struct.coordz [i1];

    v2.x = coord_vel_force_struct.coordx [i3] - coord_vel_force_struct.coordx [i2];
    v2.y = coord_vel_force_struct.coordy [i3] - coord_vel_force_struct.coordy [i2];
    v2.z = coord_vel_force_struct.coordz [i3] - coord_vel_force_struct.coordz [i2];

    v3.x = coord_vel_force_struct.coordx [i4] - coord_vel_force_struct.coordx [i3];
    v3.y = coord_vel_force_struct.coordy [i4] - coord_vel_force_struct.coordy [i3];
    v3.z = coord_vel_force_struct.coordz [i4] - coord_vel_force_struct.coordz [i3];

    _Dcoordinate m = cross_product (v1, v2);
    _Dcoordinate n = cross_product (v2, v3);

    double tmp = sqrt (v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
    double kosinus = m.x*n.x + m.y*n.y + m.z*n.z;
    double sinus = (v1.x*n.x + v1.y*n.y + v1.z*n.z) * tmp;
    return atan2 (sinus, kosinus);
}

/////////////////////////////////////////////
_Dcoordinate gen_perpendicular_vector (const _Dcoordinate &n) {
    _Dcoordinate res;
    res.x = (double) rand () / RAND_MAX;
    double p = (double) rand () / RAND_MAX;
    if (p > 0.5) res.x = -res.x;

    res.y = (double) rand () / RAND_MAX;
    p = (double) rand () / RAND_MAX;
    if (p > 0.5) res.y = -res.y;

    res.z = (double) rand () / RAND_MAX;
    p = (double) rand () / RAND_MAX;
    if (p > 0.5) res.z = -res.z;

    p = res.x*n.x + res.y*n.y + res.z*n.z;
    res.x -= p*n.x;      res.y -= p*n.y;      res.z -= p*n.z;

    p = sqrt (res.x*res.x + res.y*res.y + res.z*res.z);
    res.x /= p;          res.y /= p;          res.z /= p;
    return res;
}

////////////////////////////////////////////////////
/*void half_shift (_Dcoordinate &coord,
                 const _Dcoordinate &box) {
    if (coord.x >   0.5*box.x) coord.x -= box.x;
    if (coord.x < - 0.5*box.x) coord.x += box.x;

    if (coord.y >   0.5*box.y) coord.y -= box.y;
    if (coord.y < - 0.5*box.y) coord.y += box.y;

    if (coord.z >   0.5*box.z) coord.z -= box.z;
    if (coord.z < - 0.5*box.z) coord.z += box.z;
}*/

/////////////////////////////////////////////////////
void half_shift (double * r, double * box) {
    if      (r [0] >   0.5*box[0]) r[0] = r[0] - box[0];
    else if (r [0] < - 0.5*box[0]) r[0] = r[0] + box[0];

    if      (r [1] >   0.5*box[1]) r[1] = r[1] - box[1];
    else if (r [1] < - 0.5*box[1]) r[1] = r[1] + box[1];

    if      (r [2] >   0.5*box[2]) r[2] = r[2] - box[2];
    else if (r [2] < - 0.5*box[2]) r[2] = r[2] + box[2];
}

//////////////////////////////////////////////////////
double periodic_dist (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double *box) {
    double vec [3];
    vec [0] = x1 - x2;        vec [1] = y1 - y2;         vec[2] = z1 - z2;
    half_shift (vec, box);
    return sqrt (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

////////////////////////////////////////////
void CC_vector (int i1, int i2, double * box,
                const _coord_vel_force_struct &coord_vel_force_struct,
                double *vec) {
    vec[0] = coord_vel_force_struct.coordx [i1] - coord_vel_force_struct.coordx [i2];
    vec[1] = coord_vel_force_struct.coordy [i1] - coord_vel_force_struct.coordy [i2];
    vec[2] = coord_vel_force_struct.coordz [i1] - coord_vel_force_struct.coordz [i2];
    half_shift (vec, box);
}

////////////////////////////////////////////////
int Kdelta (int i, int j) {
    if (i != j) return (0);
    else        return (1);
}

////////////////////////////////////////////////////////////////
void find_a_name_for_this (double A, double B, double C, double H,
                           _Dcoordinate &ep1, _Dcoordinate &ep2, _Dcoordinate &ep3) {

    double a0 = sqrt (B*B - 2.0*A*B + A*A + 4.0*H*H);

    double x_max = 0.5*(B + A + a0);
    double y_max = 0.5*(B + A - a0);
    double z_max = C;

    _Dcoordinate n1, n2, n3;
    n1.x = 0.5*(B - A - a0) / H;
    n1.y = 1;
    n1.z = 0;

    double norm = sqrt (n1.x*n1.x + n1.y*n1.y + n1.z*n1.z);
    n1.x /= norm;      n1.y /= norm;      n1.z /= norm;

    n2.x = 0.5*(B - A + a0) / H;
    n2.y = 1;
    n2.z = 0;

    norm = sqrt (n2.x*n2.x + n2.y*n2.y + n2.z*n2.z);
    n2.x /= norm;      n2.y /= norm;      n2.z /= norm;

    n3.x = 0;          n3.y = 0;          n3.y = 1;

    /////////////// Sorting: x_min, y_min and z_min to keep track of ep
    double x_min, y_min;
    if (x_max < y_max) {
        x_min = x_max;
        y_min = y_max;
        ep1 = n1;
        ep2 = n2;
    } else {
        x_min = y_max;
        y_min = x_max;
        ep1 = n2;
        ep2 = n1;
    }

//    double z_min;
    if (z_max < x_min) {
//        z_min = y_min;        y_min = x_min;        x_min = z_max;
        ep3 = ep2;        ep2 = ep1;        ep1 = n3;
    } else if (z_max < y_min) {
//        z_min = y_min;        y_min = z_max;
        ep3 = ep2;        ep2 = n3;
    } else {
//        z_min = z_max;
        ep3 = n3;
    }

    if (ep1.x < 0) {
        ep1.x = -ep1.x;        ep1.y = -ep1.y;        ep1.z = -ep1.z;
    }
    if (ep2.y < 0) {
        ep2.x = -ep2.x;        ep2.y = -ep2.y;        ep2.z = -ep2.z;
    }
}

//////////////////////////////////////////////////
void compute_ep_vector (double A, double B, double C,
                        double F, double G, double H,
                        _Dcoordinate &ep1, _Dcoordinate &ep2, _Dcoordinate &ep3) {

    if (fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6)
        find_a_name_for_this (A, B, C, H, ep1, ep2, ep3);

    else if (fabs (F) < 1.0e-6 && fabs (H) < 1.0e-6)
        find_a_name_for_this (A, C, B, G, ep1, ep2, ep3);

    else if (fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6)
        find_a_name_for_this (C, B, A, F, ep1, ep2, ep3);

    else if (fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6) {
        double x_max = A, y_max = B, z_max = C;

        _Dcoordinate n1, n2, n3;
        n1.x = 1;      n1.y = 0;       n1.z = 0;
        n2.x = 0;      n2.y = 1;       n2.z = 0;
        n3.x = 0;      n3.y = 0;       n3.z = 1;

        double x_min, y_min;
        if (x_max < y_max) {
            x_min = x_max;
            y_min = y_max;
            ep1 = n1;
            ep2 = n2;
        } else {
            x_min = y_max;
            y_min = x_max;
            ep1 = n2;
            ep2 = n1;
        }

//        double z_min;
        if (z_max < x_min) {
//            z_min = y_min;            y_min = x_min;            x_min = z_max;
            ep3 = ep2;            ep2 = ep1;            ep1 = n3;
        } else if (z_max < y_min) {
//            z_min = y_min;            y_min = z_max;
            ep3 = ep2;            ep2 = n3;
        } else {
//            z_min = z_max;
            ep3 = n3;
        }

        if (ep1.x < 0) {
            ep1.x = -ep1.x;            ep1.y = -ep1.y;            ep1.z = -ep1.z;
        }

        if (ep2.y < 0) {
            ep2.x = -ep2.x;            ep2.y = -ep2.y;            ep2.z = -ep2.z;
        }

    } else {
        double a2 = - A - B - C;
        double a1 = A*B + A*C + B*C - H*H - G*G - F*F;
        double a0 = A*F*F + B*G*G + C*H*H + 2*H*G*F - A*B*C;

        double x_min = - ONE_THIRD * a2;
        double x_max = (9.0*a1*a2 - 27*a0 - 2.0*a2*a2*a2) / 54.0;

        double y_max = x_min*x_min - ONE_THIRD * a1;
        double y_min = 2*sqrt (y_max);

        double z_min = ONE_THIRD * acos (x_max * pow (y_max, - 1.5));

        a0 = ONE_THIRD * PI;

        x_max = x_min + y_min * cos (z_min);
        y_max = x_min - y_min * cos (z_min - a0);
        double z_max = x_min - y_min * cos (z_min + a0);

        if (x_max < y_max) {
            x_min = x_max;            y_min = y_max;
        } else {
            x_min = y_max;            y_min = x_max;
        }

        if (z_max < x_min) {
            z_min = y_min;            y_min = x_min;            x_min = z_max;
        } else if (z_max < y_min) {
            z_min = y_min;            y_min = z_max;
        } else z_min = z_max;

        x_max = (A - x_min)*F + G*H;
        y_max = (B - x_min)*G + F*H;
        z_max = (C - x_min)*H + G*F;

        a0 = x_max / y_max;
        a1 = x_max / z_max;
        a2 = sqrt (1.0 + a0*a0 + a1*a1);

        _Dcoordinate n1, n2, n3;
        n1.x = 1.0 / a2;
        n1.y = a0 / a2;
        n1.z = a1 / a2;

        ep1 = n1;

        x_max = (A - y_min)*F + G*H;
        y_max = (B - y_min)*G + F*H;
        z_max = (C - y_min)*H + G*F;

        a0 = y_max / x_max;
        a1 = y_max / z_max;
        a2 = sqrt (1.0 + a0*a0 + a1*a1);

        n2.x = a0 / a2;
        n2.y = 1.0 / a2;
        n2.z = a1 / a2;

        ep2 = n2;
        ep3 = n3;
    }
    ep3 = cross_product (ep1, ep2);
}

////////////////////////////////////////////////////
double comp_distance (double x1, double y1, double z1,
                      double x2, double y2, double z2) {
    double x = x1 - x2;
    double y = y1 - y2;
    double z = z1 - z2;
    return sqrt(x*x + y*y + z*z);
}

#endif
