#include <stdio.h>
#include <vector>
#include <math.h>
#include <map>
#include <utility>
#include <iostream>
#include <algorithm>

#include "common.h"

#ifndef DISPERSION_PAULI_H
#define DISPERSION_PAULI_H

class DispersionPauli:Nonbonded {

    public:
        DispersionPauli(const int num_sites,
                        const int* nuclei,
                        const double* exponents,
                        const double* coeff);
        ~DispersionPauli();

        void set_dispersion_params(double s6, double s1, double s2);
        void get_dispersion_params(double s6, double& s1, double& s2);
        void set_vdw_radii(std::map<int, double> nuclei2radiiMap);
        std::map<int, double> get_vdw_radii();

        double calc_pair_interaction(int i, int j);

    private:
        vec_i nuclei;
        vec_d exponents;
        vec_d coeff;
        vec_d vdw_radii;
        double disp_s6, disp_s1, disp_s2;
        int n_sites;

        std::vector<vec_i> exclusions;
};


#endif // DISPERSION_PAULI_H