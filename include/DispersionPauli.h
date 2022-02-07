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
        map_id get_vdw_radii_map();
        map_id get_C6_map();
        vec_d get_pauli_radii();
        vec_d get_pauli_exp();
        vec_d get_C6_coeff();
        vec_d get_vdw_radii();
        double get_pauli_energy();
        double get_disp_energy();
        int get_num_sites();
        std::vector<vec_d> get_forces();
        std::set<int> get_exclusions(const int particle1) const;

        void initialize();
        double calc_energy(const vec_d &coords);
        double calc_one_pair(const vec_d &pos, DeltaR &deltaR, int i, int j, Energies& energies);
        Energies calc_one_pair(const vec_d &pos, int i, int j);

        void create_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, int bond_cutoff);
        void create_exclusions_from_fragment(const vec_i frag_idx);
        void add_exclusion(const int i, const int j);

        void set_use_secondary_radii(bool use_radii=true);

        //  periodic boundary conditions
        void set_use_PBC(bool is_periodic);
        void set_use_PBC(bool is_periodic, const double x, const double y, const double z);
        bool get_use_PBC();

        //  pseudo-spheroid repulsion
        void set_use_two_site_repulsion(bool on_off);
        void set_two_site_distance(double vertical_dist);
        void create_repulsion_sites(double vertical_dist, const std::vector<std::pair<int, int> > &bonds);
        
    private:
        vec_i nuclei;
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

        //std::vector<Vec3> forces;
        std::vector<vec_d> forces;
        void add_force(vec_d &force, const Vec3 &dR);

        void set_all_vdw_radii();
        void set_all_C6_coeff();
        double radii_to_coeff(double radii, double exponent);
        const double kcal = 4.184/2625.5009; // 1 kcal/mol in atomic units

        bool use_secondary_radii;
        vec_d secondary_radii;
        double secondary_exp;
        map_id secondary_radii_map;
        void set_all_secondary_radii();
        Periodicity periodicity;

        //  out of plane repulsion
        std::vector<std::pair<int, int>> two_site_indicies;
        bool two_site_indicies_set, use_two_site_repulsion;
        void calc_two_site_repulsion(const vec_d &pos, DeltaR &deltaR, int i, int j, Energies& energies);
        double two_site_dist;
        Vec3 get_perp_vector(const vec_d &pos, int i, int j, int k);

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