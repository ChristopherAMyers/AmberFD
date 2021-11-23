#include <stdio.h>
#include <vector>
#include <math.h>
#include <map>
#include <utility>
#include <set>
#include <iostream>
#include <algorithm>
#include <lapacke.h>
#include <cblas.h>
#include <omp.h>

#include "common.h"
#include "Vec3.h"

#ifndef FLUC_DENS_H
#define FLUC_DENS_H

class FlucDensEnergies{
    public:
        double frz, pol, vct;
        double elec_elec, elec_nuc, nuc_nuc;
        FlucDensEnergies(){   reset();  }
        void reset(){   elec_elec = elec_nuc = nuc_nuc = frz = pol = vct = 0.0; }
        double total() { return frz + pol + vct;   }
};
        

//typedef std::vector<double> vec_d;
class FlucDens {

    public:
        FlucDens(const int n_sites, 
                   const double *frozen_charges_in, 
                   const double *nuclei_in, 
                   const double *frozen_exp, 
                   const double *dynamic_exp);
        ~FlucDens();

        void print_params(const std::string message, const std::string param_name);
        static double dot3Vec(const vec_d &coords, int i, int j);

        double calc_overlap(const vec_d &coords);
        double calc_energy(const vec_d &coords, bool calc_frz=true, bool calc_pol=true);
        void calc_one_electro(double* deltaR, int i, int j, bool calc_pol, bool calc_frz, Energies& energies);
        void initialize_calculation();
        void solve_minimization();
        void set_dampening(double coeff, double exponent);

        //  frozen - frozen exclusions
        void add_frz_frz_exclusion(int frz_i, int frz_j);
        void get_frz_frz_exclusions(const int particle1, std::vector<int> &particles2) const;
        void create_frz_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, int bond_cutoff);
        int get_num_frz_frz_exclusions() const;

        
        //  delta_rho - frozen exclusions and fragments
        void add_del_frz_exclusion(int delta_i, int frz_j);
        void add_fragment(const std::vector<int> site_idx_list);
        void change_dyn_exp(const int index, const double value);
        void change_frz_exp(const int index, const double value);
        void get_del_frz_exclusions(const int particle1, std::vector<int> &particles2) const;

        vec_d get_rho_coulomb_mat();
        vec_d get_rho_pot_vec();
        vec_d get_delta_rho();
        double get_frozen_energy();
        double get_polarization_energy();
        FlucDensEnergies get_energies();
        vec_d get_params_by_name(const std::string param_name);
        std::vector<std::vector<int>> get_constraints();

        std::vector<int> site_frag_ids;
        vec_d A_mat_save, B_vec_save;

        void set_frag_constraints(const bool constr_frags);
        double get_total_time();

    private:
        vec_d frozen_pop;
        vec_d fluc_pop;
        vec_d nuclei;
        vec_d frozen_chg;
        vec_d frozen_exp;
        vec_d dynamic_exp;
        vec_d delta_rho;
        vec_d damp_sum;
        double damp_exponent;
        double damp_coeff;
        std::vector<std::vector<int>> exclusions_del_frz;
        std::vector<std::set<int>> exclusions_frz_frz;
        size_t n_sites;
        std::map<std::string, vec_d* > param_data;
        
        std::vector<std::vector<int>> constraints;
        int n_fragments;
        bool use_frag_constraints;

        double frz_frz_overlap(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br);
        double elec_elec_energy(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br);
        double elec_nuclei_energy(const double inv_r, const double a, const double exp_ar);
        bool use_long_range_approx(double r, double a, double b);
        void create_del_exclusions_from_fragment(const std::vector<int> frag_idx);
        
        
        double dens_cutoff_pct_error = 0.02;
        double dens_cutoff_power_law;
        const double dens_cutoff_a = -1.1724, dens_cutoff_b=14.692;

        vec_d J_mat;
        vec_d pot_vec;
        vec_d dJ_dR;
        vec_d dPot_dR;
        //double total_frz_energy;
        //double total_pol_energy;

        void assign_constraints();
        std::out_of_range out_of_bounds_eror(const char *msg, const int idx1);
        std::out_of_range out_of_bounds_eror(const char *msg, const int idx1, const int idx2);
        double total_time;

        //double total_elec_elec, total_elec_nuc, total_nuc_nuc;
        
        FlucDensEnergies total_energies;
};



#endif  // FLUC_DENS_H