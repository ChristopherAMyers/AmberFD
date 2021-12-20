#include <stdio.h>
#include <vector>
#include <math.h>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <memory>

#include "common.h"
#include "FlucDens.h"
#include "DispersionPauli.h"

#ifndef AMBERFD_H
#define AMBERFD_H

/*  Small structure used to assign per-site force field parameters */
class ParticleInfo{
    public:
        int nuclei;
        double frz_chg, frz_exp, dyn_exp, pauli_exp, pauli_radii;
        ParticleInfo(int nuclei)
        {
            this->nuclei = nuclei;
            frz_chg = frz_exp = dyn_exp = pauli_exp =  pauli_radii = 0.0;
        }
        ParticleInfo(int nuclei, double frz_chg, double frz_exp, double dyn_exp, double pauli_exp, double pauli_coeff):
        nuclei(nuclei), frz_chg(frz_chg), frz_exp(frz_exp), dyn_exp(dyn_exp), pauli_exp(pauli_exp), pauli_radii(pauli_coeff)
        {}
};

class AmberFD{
    public:
        AmberFD();
        AmberFD(const int n_sites);
        ~AmberFD();

        void add_particle(ParticleInfo parameters);
        void add_fragment(const vec_i frag_idx);
        int get_num_particles()
        {   return n_sites; }
        std::vector<vec_d> get_forces();
        Energies calc_energy_forces(const vec_d &positions);
        std::shared_ptr<FlucDens> get_fluc_dens_force(bool create_if_null=false);
        std::shared_ptr<DispersionPauli> get_disp_pauli_force(bool create_if_null=false);
        std::shared_ptr<FlucDens> create_fluc_dens_force();
        std::shared_ptr<DispersionPauli> create_disp_pauli_force();

    private:
        vec_i nuclei;
        vec_d frz_chg;
        vec_d frz_exp;
        vec_d dyn_exp;
        vec_d pauli_exp;
        vec_d pauli_radii;
        int n_sites;

        std::shared_ptr<FlucDens> flucDens;
        std::shared_ptr<DispersionPauli> dispersionPauli;

        //void zero_energies();
        Energies total_energies;
        double E_total;
        std::vector<vec_d> forces;
};



#endif // AMBERFD_H