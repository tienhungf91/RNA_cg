#ifndef IO_H
#define IO_H

////////////////////
// Read maxi key
/////////////////////
void read_maxi_key (char * file_maxi_key,
                    _simu_struct  &simu_struct,
                    _topol_struct &topol_struct) {
    FILE * f1;
    f1 = fopen (file_maxi_key, "r");
    fseek (f1, 0L, SEEK_SET);

    for (int i = 1; i <= NA; i++) {
        int j;
        fscanf (f1, "%d %le %le %le %le\n", &j, topol_struct.MASS + i, topol_struct.RLJ + i, topol_struct.ELJ + i, topol_struct.QC + i);

        topol_struct.VISC [i] = simu_struct.visc_gamma * topol_struct.RLJ [i];
        topol_struct.K1   [i] = exp (-topol_struct.VISC [i] / topol_struct.MASS [i] * simu_struct.dt);
        topol_struct.K2   [i] = (1 - topol_struct.K1 [i]) / topol_struct.VISC [i];
        topol_struct.SIGMA_FORCE [i] = sqrt (2* topol_struct.VISC [i] * simu_struct.temp / simu_struct.dt);
    }
    fclose (f1);

    //// Scale down phosphate charge
    topol_struct.QC [1] *= simu_struct.Ze;

    for (int i = 1; i <= NA; i++)
        for (int j = i; j <= NA; j++) {
            int k = (i-1) * NA - (i-1)*i / 2 + j;

            // Keep all excluded interaction D0 = 3.2 A between bases
            // Check those HARD-CODE
            if (i <= NS && i >= 3 && j <= NS && j >= 3)
                topol_struct.D_LJ [k] = 3.2;
            else topol_struct.D_LJ [k] = topol_struct.RLJ [i] + topol_struct.RLJ [j];

            topol_struct.D2_LJ [k] = topol_struct.D_LJ [k] * topol_struct.D_LJ [k];

            /////////// LJ cutoff for WCA potential
            topol_struct.D2_LJ_CUTOFF [k] = topol_struct.D2_LJ [k];

            double overlap_tmp = topol_struct.D_LJ [k] - 0.2*1.5852;
            topol_struct.D2_LJ_OVERLAP [k] = overlap_tmp * overlap_tmp;

//            double min_tmp = topol_struct.D_LJ [k] - 0.4406141861*1.5852;
//            topol_struct.D2_LJ_MINIMUM [k] = min_tmp * min_tmp;

            topol_struct.E_LJ [k] = sqrt (topol_struct.ELJ [i] * topol_struct.ELJ [j]);
//            topol_struct.F_LJ [k] = 12.0 * topol_struct.E_LJ [k] / (1.5852*1.5852);
            topol_struct.F_LJ [k] = 12.0 * topol_struct.E_LJ [k];

            topol_struct.Q_C  [k] = topol_struct.QC [i] * topol_struct.QC [j];
            topol_struct.E_DH [k] = simu_struct.lB_T * topol_struct.Q_C [k];
        }
}

///////////////////////////////////////////////////////////////////////////
/*void read_pair_potential (const std::string &u_dir,
                          _simu_struct &simu_struct) {
    simu_struct.pair_potential = (double **) calloc (10, sizeof(double*));
    for (int i = 1; i <= 9; i++) {
        std::string filename;
        if      (i == 1) filename = u_dir + "/uvv.K_K";
        else if (i == 2) filename = u_dir + "/uvv.K_Cl";
        else if (i == 3) filename = u_dir + "/uvv.K_Mg";
        else if (i == 4) filename = u_dir + "/uvv.K_P";
        else if (i == 5) filename = u_dir + "/uvv.Mg_Mg";
        else if (i == 6) filename = u_dir + "/uvv.Mg_Cl";
        else if (i == 7) filename = u_dir + "/uvv.Mg_P_implicit";
        else if (i == 8) filename = u_dir + "/uvv.Cl_Cl";
        else if (i == 9) filename = u_dir + "/uvv.Cl_P";

        std::ifstream f (filename.c_str());
        if (f.is_open ()) {
            std::string line;
            int count = 0;
            while (getline(f, line)) {
                count++;
                simu_struct.pair_potential [i] = (double *) realloc (simu_struct.pair_potential [i], (count + 1) * sizeof(double));
                // Keep starting values in index 0
                if (count == 1) simu_struct.pair_potential [i][0] = atof (line.substr(0, 6).c_str());
                simu_struct.pair_potential [i][count] = atof (line.substr(7, 15).c_str());
                //printf ("%f\n", simu_struct.pair_potential [i][count]);
            }
            f.close();
        } else {
            std::cout << "Not able to read pair potential file " << filename << std::endl;
            exit (0);
        }
    }
}*/

////////////////////////////////////////////////////////////////////////
void read_pair_potential (const std::string &file,
                          _simu_struct &simu_struct) {
    std::ifstream f (file.c_str());
    if (f.is_open ()) {
        std::string line;
        int count = 0;
        while (getline(f, line)) {
            count++;
            std::stringstream stream (line);
            double r, pot;
            stream >> r >> pot;
            // Keep starting values in index 0
            if (count == 1) simu_struct.pair_potential.push_back (r);
            simu_struct.pair_potential.push_back (pot);
        }
        f.close();
    } else {
        std::cout << "Not able to read pair potential file " << file << std::endl;
        exit (0);
    }
}

//////////////////////////////////////////////////////////////
void write_rst (char * file_rst, int Natm,
                double * box,
                const _coord_vel_force_struct &coord_vel_force_struct) {
    std::ofstream f_rst (file_rst);
    if (f_rst.is_open()) {
        f_rst << "ATM " << Natm << std::endl;
        for (int i = 1; i <= Natm; i++) {
            f_rst << std::setw(16) << std::setprecision(8) << std::fixed << std::scientific << coord_vel_force_struct.coordx [i];
            f_rst << std::setw(16) << std::setprecision(8) << std::fixed << std::scientific << coord_vel_force_struct.coordy [i];
            f_rst << std::setw(16) << std::setprecision(8) << std::fixed << std::scientific << coord_vel_force_struct.coordz [i];

            f_rst << std::setw(16) << std::setprecision(8) << std::fixed << std::scientific << coord_vel_force_struct.velox [i];
            f_rst << std::setw(16) << std::setprecision(8) << std::fixed << std::scientific << coord_vel_force_struct.veloy [i];
            f_rst << std::setw(16) << std::setprecision(8) << std::fixed << std::scientific << coord_vel_force_struct.veloz [i] << std::endl;
        }

        f_rst << "BOX" << std::setw(12) << std::setprecision(4) << std::fixed << box[0];
        f_rst <<          std::setw(12) << std::setprecision(4) << std::fixed << box[1];
        f_rst <<          std::setw(12) << std::setprecision(4) << std::fixed << box[2] << std::endl;
        f_rst.close();
    }
}

///////////////////////////////////////////////////////////
void read_rst_file (char * file_rst, int Natm,
                    double * box,
                    _coord_vel_force_struct &coord_vel_force_struct) {
    std::ifstream f_rst (file_rst);
    if (f_rst.is_open()) {
        std::string line;
        int count = 0;
        int Natom_rst;
        while (getline(f_rst, line))
            if (line.find("ATM") == 0) {
                Natom_rst = atoi (line.substr(4).c_str());
                if (Natom_rst != Natm) {
                    printf ("Number of atoms not matched in restart file\n");
                    exit (1);
                }
            } else if (line.find("BOX") == 0) {
                box[0] = atof (line.substr( 6, 12).c_str());
                box[1] = atof (line.substr(18, 12).c_str());
                box[2] = atof (line.substr(30, 12).c_str());
            } else {
                count++;
                coord_vel_force_struct.coordx [count] = atof (line.substr( 1, 16).c_str());
                coord_vel_force_struct.coordy [count] = atof (line.substr(17, 16).c_str());
                coord_vel_force_struct.coordz [count] = atof (line.substr(33, 16).c_str());

                coord_vel_force_struct.velox  [count] = atof (line.substr(49, 16).c_str());
                coord_vel_force_struct.veloy  [count] = atof (line.substr(65, 16).c_str());
                coord_vel_force_struct.veloz  [count] = atof (line.substr(81, 16).c_str());
            }

        if (count != Natm) {
            printf ("Number of atoms not matched in restart file\n");
            printf ("Restart file                      %d\n", Natom_rst);
            printf ("Number of atoms read               %d\n", count);
            printf ("Number of atoms supposed to read   %d\n", Natm);
            exit (1);
        }
        f_rst.close();
    }
}

////////////////////////////////////////////////////////////
void read_progress (char * file_info,
                    long long &progress,
                    int step_traj) {
    std::ifstream f_info (file_info);
    long long tmp = 0;
    if (f_info.is_open()) {
        std::string line;
        while (getline(f_info, line))
            if (line.find(" Step") == 0)
                tmp = atoll (line.substr(5).c_str());
        f_info.close();
    }
    long int t = (long int) (tmp / step_traj);
    progress = t*step_traj;
    printf ("Restart at step   %lld\n", progress);
}

///////////////////////////////////////////////////
void write_info (char * file_info,
                 long int step,
                 double &dtime1, double &dtime2,
                // double &time_running,
                 const std::vector<int> &bead_zeroQ,
                 double &Rg,
                 const _simu_struct    &simu_struct,
                 const _crowder_struct &crowder_struct,
                 const _energy_struct  &energy_struct) {
    FILE * f_info;
    f_info = fopen (file_info, "w");
    if (f_info) {
        fprintf (f_info, " Run options\n");
        fprintf (f_info, "    -T (temp)   %10.2f   K\n", simu_struct.temp / KELVIN_TO_KT);
        if (simu_struct.pull_force) {
            fprintf (f_info, "    -f (force)  %10.4f   pN\n", simu_struct.force);
            fprintf (f_info, "    -F (f_dir)  %8.2f  %8.2f  %8.2f\n", simu_struct.f_direction.x, simu_struct.f_direction.y, simu_struct.f_direction.z);
            fprintf (f_info, "    atom pulled %8d  %8d\n", simu_struct.f_id1, simu_struct.f_id2);
        }
        fprintf (f_info, "    -M (di- or tri-valent)   %10.2e   M\n", crowder_struct.conc [0]);
        fprintf (f_info, "    -K (monovalent)          %10.2e   M\n", crowder_struct.conc [2]);
        fprintf (f_info, "          Phosphate charge     %6.3f\n", -simu_struct.Ze);
        if (bead_zeroQ.size() > 0) {
            fprintf (f_info, "          Set q=0 for those phosphates  ");
            for (int i = 0; i < bead_zeroQ.size(); i++)
                fprintf (f_info, "%6d", bead_zeroQ[i]);
            fprintf (f_info, "\n");
        }

        fprintf (f_info, "    -b (box)    %10.4f   %10.4f   %10.4f\n", simu_struct.box [0], simu_struct.box [1], simu_struct.box [2]);
        fprintf (f_info, "    -s (seed)   %d\n", simu_struct.seed);
        fprintf (f_info, "       Radius of gyration     %f  A\n", Rg);
        fprintf (f_info, "\n Step      %lld\n", simu_struct.progress);
        if (simu_struct.fix_solute == 0) {
            fprintf (f_info, "        E_BOND            %10.4f\n", energy_struct.E_BOND);
            fprintf (f_info, "        E_ANGLE           %10.4f\n", energy_struct.E_ANGLE);
            fprintf (f_info, "        E_EXCLUDED        %10.4f\n", energy_struct.E_EXCLUDED);
            fprintf (f_info, "        E_NON-NAT_HB      %10.4f\n", energy_struct.E2_HB);
            fprintf (f_info, "        E_NAT_HB          %10.4f\n", energy_struct.E3_HB);
            fprintf (f_info, "        E_NON-NAT_STACK   %10.4f\n", energy_struct.E2_STACK);
            fprintf (f_info, "        E_NAT_STACK       %10.4f\n", energy_struct.E3_STACK);
            fprintf (f_info, "        E_Q               %10.4f\n", energy_struct.E_Q);
            fprintf (f_info, "        E_potential       %10.4f\n", energy_struct.E_potential);
        }
        double step_per_day1 = step / dtime1 * 86400;
        double step_per_day2 = simu_struct.step_energy / dtime2 * 86400;
        double time_remaining = (simu_struct.nstep - step) / step_per_day1;
        //double time_remaining = time_running - dtime1;

        fprintf (f_info, "\n   Speed last %6d steps     %7.2e    steps / day\n", simu_struct.step_energy, step_per_day2);
        fprintf (f_info, "   Overall speed               %7.2e    steps / day\n", step_per_day1);
        fprintf (f_info, "   Time remaining   ");
        if (time_remaining > 1)
            fprintf (f_info, "%7.1f     days\n",  time_remaining);
        else if (time_remaining * 24 > 1)
            fprintf (f_info, "%7.1f     hours\n", time_remaining * 24);
        else fprintf (f_info, "%7.1f     mins\n", time_remaining * 1440);
        fflush (f_info);
        fclose (f_info);
    }
}

#endif
