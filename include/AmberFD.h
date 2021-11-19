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

class ParticleInfo{
    public:
        int nuclei;
        double frz_chg, frz_exp, dyn_exp, pauli_exp, pauli_coeff;
        ParticleInfo(int nuclei)
        {
            this->nuclei = nuclei;
            frz_chg = frz_exp = dyn_exp = pauli_exp =  pauli_coeff = 0.0;
        }
        ParticleInfo(int nuclei, double frz_chg, double frz_exp, double dyn_exp, double pauli_exp, double pauli_coeff):
        nuclei(nuclei), frz_chg(frz_chg), frz_exp(frz_exp), dyn_exp(dyn_exp), pauli_exp(pauli_exp), pauli_coeff(pauli_coeff)
        {}
};

class AmberFD{
    public:
        AmberFD();
        AmberFD(const int n_sites);
        ~AmberFD();

        void add_particle(ParticleInfo parameters);
        std::shared_ptr<FlucDens> create_fluc_dens_force();
        std::shared_ptr<DispersionPauli> create_disp_pauli_force();

    private:
        vec_i nuclei;
        vec_d frz_chg;
        vec_d frz_exp;
        vec_d dyn_exp;
        vec_d pauli_exp;
        vec_d pauli_coeff;
        int n_sites;

        std::shared_ptr<FlucDens> flucDens;
        std::shared_ptr<DispersionPauli> dispersionPauli;
};



#endif // AMBERFD_H