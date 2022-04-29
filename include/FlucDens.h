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
        void calc_one_electro(DeltaR &deltaR, int i, int j, bool calc_pol, bool calc_frz, Energies& energies, std::vector<Vec3> &forces, int thread_num=0);
        Energies calc_one_frozen(const vec_d &coords, int i, int j);
        void initialize_calculation();
        void solve_minimization(std::vector<Vec3> &forces);
        //void set_dampening(double coeff, double exponent);
        void set_dampening(double coeff, double exponent, double pol_wall=1.2);

        //  frozen - frozen exclusions
        void add_frz_frz_exclusion(int frz_i, int frz_j);
        std::set<int>  get_frz_frz_exclusions(const int particle1) const;
        void create_frz_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, int bond_cutoff);
        int get_num_frz_frz_exclusions() const;

        
        //  delta_rho - frozen exclusions and fragments
        void add_del_frz_exclusion(int delta_i, int frz_j);
        std::set<int> get_del_frz_exclusions(const int particle1) const;
        void add_fragment(const std::vector<int> site_idx_list);
        int get_num_fragments();
        std::vector<vec_i> get_fragments();
        
        void set_dyn_exp(const int index, const double value);
        void set_dyn_exp(vec_d exponents);
        void set_frz_exp(const int index, const double value);
        void set_ct_coeff(const double coeff);
        void set_site_params(const int index, const double frz_chg, const double frz_exp, const double dyn_exp);
        void set_additional_hardness(const int index, const double value);
        void set_additional_hardness(vec_d values);

        int get_num_constraints();
        vec_d get_rho_coulomb_mat();
        vec_d get_rho_pot_vec();
        vec_d get_delta_rho();
        void get_dampening(double &coeff, double &exponent);
        double get_ct_coeff();
        double get_frozen_energy();
        double get_polarization_energy();
        double get_ct_energy();
        FlucDensEnergies get_energies();
        std::vector<std::string> get_param_names();
        vec_d get_params_by_name(const std::string param_name);
        std::vector<std::vector<int>> get_constraints();
        void get_site_params(const int index, double &frz_chg, double &frz_exp, double &dyn_exp);

        std::vector<int> site_frag_ids;
        std::vector<std::vector<int> > site_frag_ids_partitioned;
        vec_d A_mat_save, B_vec_save;

        void set_frag_constraints(const bool constr_frags);
        void set_calc_forces(bool calculate_forces);
        double get_total_time();
        std::vector<vec_d> get_forces();

        //  single energy terms
        double frz_frz_overlap(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br);
        double elec_elec_energy(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br, double &dEdR);
        double elec_nuclei_energy(const double inv_r, const double a, const double exp_ar, double &dEdR);


        //  periodic boundary conditions
        void set_use_PBC(const bool is_periodic);
        void set_use_PBC(const bool is_periodic, const double x, const double y, const double z);
        bool get_use_PBC();

        //  cutoff info
        void set_cutoff_distance(double distance_in_nm);
        double get_cutoff_distance();
        void set_use_cutoff(bool useCutoff);
        bool get_use_cutoff();
        void set_use_SR_cutoff(bool);
        bool get_use_SR_cutoff();

        //  Density printing
        enum DensityType {
            All = 0, Frozen=1, Delta=2, Nuclei=3
        };
        vec_d calc_density(const vec_d &points, const vec_d &pos, DensityType density_type);

        //  External fields and dipoles
        bool has_ext_field;
        void set_external_field(double field_x, double field_y, double field_z);
        void apply_field_to_system(const vec_d &coords);
        std::vector<Vec3> get_dipoles(const vec_d &coords);
        vec_d get_dipole(const vec_d &coords, DensityType density_dype);
        double calc_frz_ext_field_energy(const vec_d &coords, std::vector<Vec3> &forces);




    private:
        
        vec_d fluc_pop;
        vec_d nuclei;
        vec_d frozen_pop;
        vec_d frozen_chg;
        vec_d frozen_exp;
        vec_d dynamic_exp;
        vec_d delta_rho;
        vec_d hardness;
        //vec_d damp_sum;

        double damp_exponent;
        double damp_coeff;
        double pol_wall_coeff;
        double pol_wall_exponent;
        double ct_coeff;
        std::vector<std::set<int>> exclusions_del_frz;
        std::vector<std::set<int>> exclusions_frz_frz;
        size_t n_sites;
        std::map<std::string, vec_d* > param_data;
        std::vector<std::vector<int>> constraints;
        int n_fragments;
        bool use_frag_constraints;
        bool remove_core;
        

        bool use_SR_approx(double r, double a, double b);
        void create_del_exclusions_from_fragment(const std::vector<int> frag_idx);

        vec_d J_mat;
        vec_d pot_vec;
        vec_d potential_mat;

        std::vector<Vec3> dJ_dPos;
        std::vector<Vec3> dDamp_dPos;
        std::vector<Vec3> dPot_dPos, dPot_dPos_trans;
        std::vector<Vec3> frozen_forces;
        std::vector<Vec3> total_forces;
        std::vector<vec_d> thread_dampening;

        void assign_constraints();
        std::out_of_range out_of_bounds_eror(const char *msg, const int idx1);
        std::out_of_range out_of_bounds_eror(const char *msg, const int idx1, const int idx2);
        double total_time;
        bool calc_forces;
        
        FlucDensEnergies total_energies;
        Periodicity periodicity;

        double cutoff_distance;
        bool use_cutoff;
        bool use_SR_cutoff;
        double SR_cutoff_pct_error = 0.10;
        const double SR_cutoff_a = -1.1724, SR_cutoff_b=14.692;
        double SR_cutoff_coeff;

        //  external fields
        
        Vec3 ext_field;
        vec_d ext_field_potential;
};



#endif  // FLUC_DENS_H