#include <stdio.h>
#include <vector>
#include <math.h>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <algorithm>

#include "common.h"

#ifndef DISPERSION_PAULI_H
#define DISPERSION_PAULI_H

class DispersionPauli {

    public:
        DispersionPauli(const int num_sites,
                        const int* nuclei,
                        const double* exponents,
                        const double* radii);
        ~DispersionPauli();

        void set_dispersion_params(double s6, double a1, double a2);
        void set_vdw_radii(map_id nuclei2radiiMap);
        void set_C6_map(map_id nucleiToC6Map);
        void set_pauli_radii(vec_d radii_list);
        void set_pauli_radii(int index, double radii);
        void set_pauli_exp(vec_d exp_list);
        void set_pauli_exp(int index, double exponent);

        void get_dispersion_params(double &s6, double &a1, double &a2);
        map_id get_vdw_radii();
        map_id get_C6_map();
        vec_d get_pauli_radii();
        vec_d get_pauli_exp();
        vec_d get_C6_coeff();
        double get_pauli_energy();
        double get_disp_energy();
        int get_num_sites();

        double calc_energy(const vec_d &coords);
        double calc_one_pair(DeltaR &deltaR, int i, int j, Energies& energies);

        void create_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, int bond_cutoff);
        void create_exclusions_from_fragment(const vec_i frag_idx);
        void add_exclusion(const int i, const int j);

        //class Info;
        vec_i nuclei;
    private:
        
        vec_d pauli_exponents;
        vec_d pauli_coeff;
        vec_d pauli_radii;
        vec_d vdw_radii;
        vec_d C6_coeff;
        double disp_s6, disp_a1, disp_a2;
        int n_sites;

        map_id C6_map;
        map_id vdw_radii_map;
        std::vector<std::set<int>> exclusions;

        double total_disp_energy;
        double total_pauli_energy;

        void set_all_vdw_radii();
        void set_all_C6_coeff();
        double radii_to_coeff(double radii, double exponent);
        const double kcal = 4.184/2625.5009; // 1 kcal/mol in atomic units

        
};

class Info
{
    public:
        double one, two, three;
        Info(){
            one = two = three = 0.0;
        }
};

#endif // DISPERSION_PAULI_H