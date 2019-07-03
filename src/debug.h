#ifndef DEBUG_H
#define DEBUG_H

//////////////////////////////////////////////////////////////////
void debug_force (const std::vector<_Res_cg> &res,
                  _simu_struct            &simu_struct,
                  _topol_struct           &topol_struct,
                  _hydbond_struct         &hydbond_struct,
                  _stack_struct           &stack_struct,
                  _coord_vel_force_struct &coord_vel_force_struct,
                  _energy_struct          &energy_struct) {
                  //const _EW_struct        &EW_struct) {

    ///// For debugging purpose and not affecting running performance
    unsigned int myseed = 1;
    FILE * f_force;
    f_force = fopen ("force.dat", "w");

    fprintf (f_force, "##### total force #####\n");
    for (int i = 1; i <= topol_struct.Natm; i++) {
        double f = coord_vel_force_struct.forcex [i]*coord_vel_force_struct.forcex [i] + coord_vel_force_struct.forcey [i]*coord_vel_force_struct.forcey [i] + coord_vel_force_struct.forcez [i]* coord_vel_force_struct.forcez [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.forcex [i], coord_vel_force_struct.forcey [i], coord_vel_force_struct.forcez [i], sqrt (f));
    }

/*    fprintf (f_force, "##### LR force #####\n");
    for (int i = 1; i <= topol_struct.Natm; i++) {
        double f = coord_vel_force_struct.flx [i]*coord_vel_force_struct.flx [i] + coord_vel_force_struct.fly [i]*coord_vel_force_struct.fly [i] + coord_vel_force_struct.flz [i]*coord_vel_force_struct.flz [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.flx [i], coord_vel_force_struct.fly [i], coord_vel_force_struct.flz [i], sqrt (f));
    }

    for (int i = 1; i <= topol_struct.Natm; i++) {
        coord_vel_force_struct.forcex [i] = 0;
        coord_vel_force_struct.forcey [i] = 0;
        coord_vel_force_struct.forcez [i] = 0;
    }

    SR_interactions_lists (simu_struct, topol_struct, coord_vel_force_struct, EW_struct);
    fprintf (f_force, "##### SR force #####\n");
    for (int i = 1; i <= topol_struct.Natm; i++) {
        double f = coord_vel_force_struct.forcex [i]*coord_vel_force_struct.forcex [i] + coord_vel_force_struct.forcey [i]*coord_vel_force_struct.forcey [i] + coord_vel_force_struct.forcez [i]* coord_vel_force_struct.forcez [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.forcex [i], coord_vel_force_struct.forcey [i], coord_vel_force_struct.forcez [i], sqrt (f));
    }*/

    for (int i = 1; i <= topol_struct.Natm; i++) {
        coord_vel_force_struct.forcex [i] = 0;
        coord_vel_force_struct.forcey [i] = 0;
        coord_vel_force_struct.forcez [i] = 0;
    }

    non_bonded_interaction (simu_struct, topol_struct, coord_vel_force_struct, energy_struct);
    fprintf (f_force, "##### non_bonded force #####\n");
    for (int i = 1; i <= topol_struct.Natm; i++) {
/*        double fx = coord_vel_force_struct.forcex [i] + coord_vel_force_struct.flx [i];
        double fy = coord_vel_force_struct.forcey [i] + coord_vel_force_struct.fly [i];
        double fz = coord_vel_force_struct.forcez [i] + coord_vel_force_struct.flz [i];
        double f = sqrt (fx*fx + fy*fy + fz*fz);*/

        double f = coord_vel_force_struct.forcex [i]*coord_vel_force_struct.forcex [i] + coord_vel_force_struct.forcey [i]*coord_vel_force_struct.forcey [i] + coord_vel_force_struct.forcez [i]* coord_vel_force_struct.forcez [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.forcex [i], coord_vel_force_struct.forcey [i], coord_vel_force_struct.forcez [i], sqrt (f));
    }

    for (int i = 1; i <= topol_struct.Natm; i++) {
        coord_vel_force_struct.forcex [i] = 0;
        coord_vel_force_struct.forcey [i] = 0;
        coord_vel_force_struct.forcez [i] = 0;
    }

    polymer_constraints (simu_struct.box, res, topol_struct, coord_vel_force_struct, energy_struct);
    fprintf (f_force, "##### polymer force #####\n");
    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        double f = coord_vel_force_struct.forcex [i]*coord_vel_force_struct.forcex [i] + coord_vel_force_struct.forcey [i]*coord_vel_force_struct.forcey [i] + coord_vel_force_struct.forcez [i]* coord_vel_force_struct.forcez [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.forcex [i], coord_vel_force_struct.forcey [i], coord_vel_force_struct.forcez [i], sqrt (f));
    }

    for (int i = 1; i <= topol_struct.Natm; i++) {
        coord_vel_force_struct.forcex [i] = 0;
        coord_vel_force_struct.forcey [i] = 0;
        coord_vel_force_struct.forcez [i] = 0;
    }

    stacking_interactions (topol_struct, coord_vel_force_struct, stack_struct, energy_struct);
    fprintf (f_force, "\n##### stack force #####\n");
    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        double f = coord_vel_force_struct.forcex [i]*coord_vel_force_struct.forcex [i] + coord_vel_force_struct.forcey [i]*coord_vel_force_struct.forcey [i] + coord_vel_force_struct.forcez [i]* coord_vel_force_struct.forcez [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.forcex [i], coord_vel_force_struct.forcey [i], coord_vel_force_struct.forcez [i], sqrt (f));
    }

    for (int i = 1; i <= topol_struct.Natm; i++) {
        coord_vel_force_struct.forcex [i] = 0;
        coord_vel_force_struct.forcey [i] = 0;
        coord_vel_force_struct.forcez [i] = 0;
    }

    hydrogen_bonds_interactions (coord_vel_force_struct, hydbond_struct, simu_struct.temp, energy_struct, &myseed);
    fprintf (f_force, "\n##### hydbond force #####\n");
    for (int i = 1; i <= topol_struct.Natm_biomol; i++) {
        double f = coord_vel_force_struct.forcex [i]*coord_vel_force_struct.forcex [i] + coord_vel_force_struct.forcey [i]*coord_vel_force_struct.forcey [i] + coord_vel_force_struct.forcez [i]* coord_vel_force_struct.forcez [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.forcex [i], coord_vel_force_struct.forcey [i], coord_vel_force_struct.forcez [i], sqrt (f));
    }

    fclose (f_force);
}

//////////////////////////////////////////////////////////////////
void debug_force_fixed (int Natm,
                        _coord_vel_force_struct &coord_vel_force_struct) {

    ///// For debugging purpose and not affecting running performance
    FILE * f_force;
    f_force = fopen ("force.dat", "w");

    fprintf (f_force, "##### total force #####\n");
    for (int i = 1; i <= Natm; i++) {
        double f = coord_vel_force_struct.forcex [i]*coord_vel_force_struct.forcex [i] + coord_vel_force_struct.forcey [i]*coord_vel_force_struct.forcey [i] + coord_vel_force_struct.forcez [i]* coord_vel_force_struct.forcez [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.forcex [i], coord_vel_force_struct.forcey [i], coord_vel_force_struct.forcez [i], sqrt (f));
    }

/*    fprintf (f_force, "##### LR force #####\n");
    for (int i = 1; i <= Natm; i++) {
        double f = coord_vel_force_struct.flx [i]*coord_vel_force_struct.flx [i] + coord_vel_force_struct.fly [i]*coord_vel_force_struct.fly [i] + coord_vel_force_struct.flz [i]*coord_vel_force_struct.flz [i];
        fprintf (f_force, "%e   %e   %e   %e\n", coord_vel_force_struct.flx [i], coord_vel_force_struct.fly [i], coord_vel_force_struct.flz [i], sqrt (f));
    }

    fprintf (f_force, "##### SR force #####\n");
    for (int i = 1; i <= Natm; i++) {
        double fx = coord_vel_force_struct.forcex [i] - coord_vel_force_struct.flx [i];
        double fy = coord_vel_force_struct.forcey [i] - coord_vel_force_struct.fly [i];
        double fz = coord_vel_force_struct.forcez [i] - coord_vel_force_struct.flz [i];
        double f = sqrt (fx*fx + fy*fy + fz*fz);
        fprintf (f_force, "%e   %e   %e   %e\n", fx, fy, fz, f);
    }*/

    fclose (f_force);
}

#endif
