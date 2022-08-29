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
#include "timer.h"
#include "LinearSolvers.h"

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
        

/**
 * This class impliments the fluctuating monopole density interaction. 
 * Similar to a fluctuating charge model, this force uses atom centered
 * monopole densities described by slater functions for the electrons and
 * point-particles for the nuclei. As such, FlucDens includes effects of
 * electron overlap between atomic sites. 
 * 
 * To use this force, create a FlucDens force object by specifying the number of 
 * charge sites and the parameters for all sites. These parameters can be altered
 * after the creation of the force object. Next, each site must be assigned a
 * molecular "fragment" Each fragment, defined by it's site indicies, are groups
 * of charge sites whose dynamic desity coefficients will sum to zero when the electrostatic
 * energy is minized. As a result, charge transfer is not allowed between fragments, but
 * each fragment will however reduce to it's frozen electron populations when 
 * complety separated from all other framgents. Each fragment is also used to
 * define which frozen density - dynamic density interactions are included. Frozen
 * densities on each site do not interact with the dynamic densities on each fragment,
 * allowing for this limmiting behavior to occur.
 * 
 * Next, set the model's pther settings by calling their respective setter function.
 * This includes the polarization dampening strengths, charge-transfer coefficients, and
 * additional hardness values to be applied to the Coulomb matrix. Polarization dampening 
 * is applied between a polarized density site and it's interacing frozen density site and
 * decays exponentially as a function of their distance. The strength and rate of this 
 * polarization can be adjusted, as well as whether this is applied to the linear  
 * (default) or the quadratic matrix terms (old method). Charge transfer is not explicitly
 * allowed in this force, and instead the energy associated with this interation is estimated
 * by a scalar multiple of the total polarization energy.
 * 
 * Additional hardness values can be applied to the Coulomb matrix. These values are not overrides, 
 * but instead add a constant value to the self-interaction between two densities and hinder 
 * (or enhance) the ability of that site to gain or loose electrons. These values can be 
 * either positive or negative, with positive values increasing the hardness. Care must be
 * taken that the system remains positive-definite if negative values are used.
 * 
 * Aconstant electric field can also be applied to each system. Simply call 
 * set_external_field() with the desired field components to apply this field. Periodic 
 * boundary conditions are also possible, and can be set by calling set_use_PBC(). If PBC
 * is used, then a long-range cutoff distance must be used (the default is to use half the
 * minimum periodic box lengths). Turning this off and using PBC will have undefined behavior.
 * 
 * This module also computes all interaction energies, including minimizing the total
 * energy w.r.t. the dynamic density coefficients. Do run this computation, first call
 * calc_energy() to compute all frozen-frozen energy components and to create the
 * polarization matrix components. Then, call solve_minimization() to perform the actual
 * matrix algebra associated with the minimization procedures. Energy components can then be
 * returned calling get_energies() or any of the particular energy getter functions.
 * Dipoles (including the polarized density) can also then be called with the updated
 * electron populations.
 * 
 * 
 * All values in this module use atomic units unless otherwise specified.
 * 
 * @brief To construct a new FlucDens force object, define the object with
 * 
 * @param n_sites number of particles (including cirtual sites)
 * @param frozen_charges_in the frozen charges of all sites
 * @param nuclei_in the nuclei of each site. Core electrons will be remove automatically 
 * @param frozen_exp the Slater exponents of all frozen densities
 * @param dynamic_exp the Slater exponents of all dynamic densities
 * 
 */
class FlucDens {
    public:
        /**
         * @brief Construct a new FlucDens force object
         * 
         * @param n_sites number of particles (including cirtual sites)
         * @param frozen_charges_in array of frozen charges of all sites
         * @param nuclei_in array of nuclei of each site. Core electrons will be remove automatically 
         * @param frozen_exp array of Slater exponents of all frozen densities
         * @param dynamic_exp array of Slater exponents of all dynamic densities
         */
        FlucDens(const int n_sites, 
                   const double *frozen_charges_in, 
                   const double *nuclei_in, 
                   const double *frozen_exp, 
                   const double *dynamic_exp);
        ~FlucDens();

        /**
         * @brief enumeration of the types of polarization dampening to be used 
         * 
         */
        enum DampType {
            /**
             * @brief Dampening will be applied to the linear terms
             */
            Linear = 1, 
            /**
             * @brief Dampening will be applied to the quadratic terms (old method)
             */
            Quadratic = 2
        };
        /**
         * @brief enumeration of the types of densities to be computed
         * 
         */
        enum DensityType {
            /**
             * @brief The total density of the system (nuclei, frozen, and dynamic)
             */
            All = 0,
            /**
             * @brief Frozen electron density only
             */
            Frozen=1,
            /**
             * @brief Dynamic electron density only
             */
            Delta=2,
            /**
             * @brief Nuclei charge
             */
            Nuclei=3
        };
        /**
         * @brief Compute the electron density at the specified coordinates
         * 
         * @param points A 3K array of K number of points to evaluate the density at
         * @param pos The 3N array of positions for the charge sites of the system
         * @param density_type The type of density to compute
         * @return vec_d A K-length array of density evaluations
         */
        vec_d calc_density(const vec_d &points, const vec_d &pos, DensityType density_type);
        //****************************************************************
        //                   indicidual site parameters 
        //****************************************************************
        /**
         * @brief Change the parameters for one charge site
         * 
         * @param index The index of the site to adjust the parameters
         * @param frz_chg The frozen charge to set
         * @param frz_exp The Slater exponent of the frozen density
         * @param dyn_exp The Slater exponent of the dynamic density
         */
        void set_site_params(const int index, const double frz_chg, const double frz_exp, const double dyn_exp);
        /**
         * @brief Get the parameters of an atomic site
         * 
         * @param index The index of the site to adjust the parameters
         * @param frz_chg The frozen charge of the site
         * @param frz_exp The Slater exponent of the frozen density
         * @param dyn_exp The Slater exponent of the dynamic density
         */
        void get_site_params(const int index, double &frz_chg, double &frz_exp, double &dyn_exp);
        /**
         * @brief Set the Slater exponent of a particular dynamic density site
         * 
         * @param index The index of the site
         * @param value The value of the Slater exponent
         */
        void set_dyn_exp(const int index, const double value);
        /**
         * @brief Set the Slater exponents of all dynamic density sites at once
         * 
         * @param exponents An array of all exponents in order of their indicies
         */
        void set_dyn_exp(vec_d exponents);
        /**
         * @brief Set the Slater exponent of a particular frozen density site
         * 
         * @param index The index of the site
         * @param value The value of the Slater exponent
         */
        void set_frz_exp(const int index, const double value);
        /**
         * @brief Set the charge-transfer coefficient. Charge transfer is extimated as
         * a multiple of the polarization energy between two inducable sites. This does
         * not apply to the polarization due to external fields
         * 
         * @param coeff The scalar multiple of the polarization energy used for the charge transfer energy
         */
        void set_ct_coeff(const double coeff);
        /**
         * @brief Get the charge-transfer coefficient
         * 
         * @return double The coefficient
         */
        double get_ct_coeff();
        /**
         * @brief Set the polarization dampening parameters
         * 
         * @param coeff The strength of the dampening  
         * @param exponent The exponential decay of the dampening
         * @param damp Tye type of dampening, enumerated by DampType
         */
        void set_dampening(double coeff, double exponent, DampType damp = DampType::Linear);
        /**
         * @brief Return the polarization dampening parameters
         * 
         * @param coeff The strength of the dampening   
         * @param exponent The exponential decay of the dampening
         */
        void get_dampening(double &coeff, double &exponent);
        /**
         * @brief Add additional atomic hardness to the diagonal J-matrix elements. This is added
         * to the value already computed by the self interaction in the matrix, not an override. 
         * 
         * @param index The particle to add the hardness to.
         * @param value The value of the additional hardness
         */
        void set_additional_hardness(const int index, const double value);
        /**
         * @brief Return the additional hardness values of all sites
         * 
         * @param values An array of the hardness values
         */
        void set_additional_hardness(vec_d values);


        //****************************************************************
        //          delta_rho - frozen exclusions and fragments 
        //****************************************************************
        /**
         * @brief Create a new fragment whoes dynamic density whould sum to zero.
         * Every site must be assigned to a fragment, else an exception will be thrown.
         * 
         * @param site_idx_list A list of site indicies that make up this fragment
         */
        void add_fragment(const std::vector<int> site_idx_list);
        /**
         * @brief Return the partitioned fragments and their site indicies for each
         * 
         * @return An array of fragment indicies. Each element correpsonds to a fragment
         * and contains a list of indicies that make up that fragment
         */
        std::vector<vec_i> get_fragments();
        /**
         * @brief Get the num fragments
         * 
         * @return The number of fragments
         */
        int get_num_fragments();
        /**
         * @brief Exclude the frozen electrostatic interactions between a dynamic and frozen site
         * 
         * @param delta_i The index of the dynamic density site
         * @param frz_j The index of the frozen density site
         */
        void add_del_frz_exclusion(int delta_i, int frz_j);
        /**
         * @brief Get the indicies of all frozen sites that do not interact with a single dynamic density
         * 
         * @param particle1 The index of the dynamic density 
         * @return std::set<int> The indicies of the frozen sites that do not interact with this
         * dynamic density
         */
        std::set<int> get_del_frz_exclusions(const int particle1) const;
        /**
         * @brief Choose whether constraints are used on all fragments or not. NOT FULLY IMPLIMENTED YET!
         * 
         * @param constr_frags If true, fragment constraints will be used. If false, then only
         * the total dynamic density of the entire system will be summed to zero, not per fragment.
         */
        void set_frag_constraints(const bool constr_frags);
        /**
         * @brief Get the number of total constrains applied to the energy minimization
         * 
         * @return The total number of constraints
         */
        int get_num_constraints();
        /**
         * @brief Get constraint matrix applied to the dyanmic densities
         * 
         * @return std::vector<std::vector<int>> 
         */
        std::vector<std::vector<int>> get_constraints();
        
        
        //****************************************************************
        //                frozen - frozen exclusions 
        //****************************************************************
        /**
         * @brief Exclude the interaction of two frozen densities
         * 
         * @param frz_i The index of the first frozen site
         * @param frz_j The index of the second frozen site
         */
        void add_frz_frz_exclusion(int frz_i, int frz_j);
        /**
         * @brief Get the all sites that are excluded from interacting from a particular site
         * 
         * @param particle1 The frozen site whoes exclusions are being inquired
         * @return std::set<int> All indicies that are excluded
         */
        std::set<int>  get_frz_frz_exclusions(const int particle1) const;
        /**
         * @brief Get the number of frz-frz exclusions 
         * 
         * @return int The total number of exclusions in the system
         */
        int get_num_frz_frz_exclusions() const;
        /**
         * @brief Create frozen-frozen exclusions based on each site's bonds.
         * 
         * @param bonds An vector of pair of indicies for each bond  
         * @param bond_cutoff All indicies that are this number of bonds away will be excluded
         */
        void create_frz_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, int bond_cutoff);
        

        //****************************************************************
        //                  External fields and dipoles 
        //****************************************************************
        /**
         * @brief Apply a constant electric field to the entire system
         * 
         * @param field_x The x-direction Cartesian component of the field
         * @param field_y The y-direction Cartesian component of the field
         * @param field_z The z-direction Cartesian component of the field
         */
        void set_external_field(double field_x, double field_y, double field_z);
        /**
         * @brief Apply the field to the system. This must be called before any minimization
         * of the total energy is performed.
         * 
         * @param coords 1D array of positions of the all sites
         */
        void apply_field_to_system(const vec_d &coords);
        /**
         * @brief Return the electric dipoles of the system. This includes nuclear, frozen electrons,
         * dynamic electrons, and total dipoles
         * 
         * @param coords 1D array of positions of the all sites
         * @return std::vector<Vec3> The dipoles that are computed. Indicies are determind by 
         * DensityType enumeration
         */
        std::vector<Vec3> get_dipoles(const vec_d &coords);
        /**
         * @brief Compute the dipole for a particule type of charge
         * 
         * @param coords 1D array of positions of the all sites
         * @param density_dype DensityType: The type of charges to be used in the computation
         * @return The dipole components
         */
        vec_d get_dipole(const vec_d &coords, DensityType density_dype);
        /**
         * @brief Compute the energy and forces of the frozen charges (nuclei and frozen electrons)
         * with the external electric field
         * 
         * @param coords 1D array of positions of the all sites
         * @param forces Array of forces to be updated
         * @return The value of the interaction energy 
         */
        double calc_frz_ext_field_energy(const vec_d &coords, std::vector<Vec3> &forces);


        //****************************************************************
        //                  Periodic boundary conditions
        //****************************************************************
        /**
         * @brief Turn on or off periodic boundary conditions
         * 
         * @param is_periodic If true, use PBC, if false, do not use them
         */
        void set_use_PBC(const bool is_periodic);
        /**
         * @brief Set the periodic boundaries
         * 
         * @param is_periodic If true, use PBC, if false, do not use them
         * @param x The x-direction of the periodic box
         * @param y The y-direction of the periodic box
         * @param z The z-direction of the periodic box
         */
        void set_use_PBC(const bool is_periodic, const double x, const double y, const double z);
        /**
         * @brief Inquire if PBC is used or not
         * 
         * @return true PBC is being used
         * @return false PBC is not being used
         */
        bool get_use_PBC();


        //****************************************************************
        //                      Cutoff distances
        //****************************************************************
        /**
         * @brief Set whether or not to use a nonbonded cutoff. Care must be taked to ensure
         * that the polarization J-matrix remains is positive definite.
         * 
         * @param useCutoff If true, enable the use of a cutoff distance between all interactions.
         * If false, do not use any cutoff distance.
         */
        void set_use_cutoff(bool useCutoff);
        /**
         * @brief Inquire of nonbonded cutoffs are being used
         * 
         * @return true cutoffs are being used
         * @return false cutoffs are not being used
         */
        bool get_use_cutoff();
        /**
         * @brief Set the cutoff distance (in bohr)
         * 
         * @param distance_in_nm The cutoff distance (in bohr)
         */
        void set_cutoff_distance(double distance_in_nm);
        /**
         * @brief Get the cutoff distance (in bohr)
         * 
         * @return distance_in_nm The cutoff distance (in bohr)
         */
        double get_cutoff_distance();
        /**
         * @brief Set the use SR cutoff. If enabled, electron overlap will be ignored
         * beyond a particular distance.
         * 
         * @param use_cutoff If true, enable the use this apprximation
         */
        void set_use_SR_cutoff(bool use_cutoff);
        /**
         * @brief Inquire of SR cutoffs are being used
         * 
         * @return true SR cutoffs are being used
         * @return false SR cutoffs are not being used
         */
        bool get_use_SR_cutoff();


        //****************************************************************
        //                      Energies and Forces
        //****************************************************************
        /**
         * @brief Calculate the total electron overlap
         * 
         * @param coords the positions of all sites
         * @param density_dype Frozen, Delta, or Total population type 
         * @return the total overlap in e^2 /bohr^3
         */
        double calc_overlap(const vec_d &coords, DensityType density_type, const vec_i &indices = vec_i());
        /**
         * @brief Calculate the energy and forces of the system.
         * 
         * @param coords A 1D array of all the positions (in bohr)
         * @param calc_frz Whether or not to compute frozen electrostatic energies
         * @param calc_pol Whether or not to compute polarization
         * @return The total frozen energy computed. Must cal solve_minimization() in 
         * order to get polarization energy.
         */
        double calc_energy(const vec_d &coords, bool calc_frz=true, bool calc_pol=true);
        /**
         * @brief Compute one pair of interactions between two atomic sites
         * 
         * @param deltaR The distance object between the two sites
         * @param i The index of the first site
         * @param j The index of the second site
         * @param calc_pol Whether or not to include polarization
         * @param calc_frz Whether or not to include frozen electrostatic energies
         * @param energies Returns the energies computed
         * @param forces Returns the forces between both sites
         * @param thread_num the index of the OpenMP thread to compute the energies and forces on
         */
        void calc_one_electro(DeltaR &deltaR, int i, int j, bool calc_pol, bool calc_frz, Energies& energies, std::vector<Vec3> &forces, int thread_num=0);
        /**
         * @brief Compute one pair of frozen interactions between two atomic sites (does not include polarization)
         * 
         * @param coords 1D array of all the positions (in bohr)
         * @param i The index of the first site
         * @param j The index of the second site
         * @return Energies energies of the computed pair
         */
        Energies calc_one_frozen(const vec_d &coords, int i, int j);
        /**
         * @brief initizlize arrays before an energy computation is performed. This is called
         * automatically by calc_energy(), but must be called manually to clear all energy and
         * force arrays if calling calc_one_electro() or calc_one_frozen().
         */
        void initialize_calculation();
        /**
         * @brief Minimize the total energy of the system w.r.t. the fluctuating densities. This
         * function solve the matrix equations that minimizes the quadratic polarization energy,
         * under the constraint that all fragment dynamic densities sum to zero.
         * 
         * @param forces array of forces to be updated based on the polarization interactions
         */
        void solve_minimization(std::vector<Vec3> &forces);
        /**
         * @brief Get the energy associated with frozen-frozen dsitity interactions
         * 
         * @return double frozen-frozen energy
         */
        double get_frozen_energy();
        /**
         * @brief Get the polarization energy
         * 
         * @return double polarization energy
         */
        double get_polarization_energy();
        /**
         * @brief Get the estimated CT energy as a multiple of the polarization energy
         * 
         * @return double 
         */
        double get_ct_energy();
        /**
         * @brief Return all energy components computed
         * 
         * @return FlucDensEnergies Object that constains the energy components
         */
        FlucDensEnergies get_energies(); 
        /**
         * @brief Get the forces exerted on all particles after an energy computation is performed
         * 
         * @return std::vector<vec_d> Nx3 Array of forces on all sites
         */
        std::vector<vec_d> get_forces();
        /**
         * @brief Overlap between frozen-frozen dsnsity
         * 
         * @param inv_r inverse distance between the two sites
         * @param a frozen exponent of the first site
         * @param b frozen exponent of the second site
         * @param exp_ar pre-computed exp(-a*r) using the first exponent
         * @param exp_br pre-computed exp(-a*r) using the secpnd exponent
         * @return double the overlap in e^2/bohr^3
         */
        double frz_frz_overlap(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br);
        /**
         * @brief Compute the density-density interaction, either frozen or dynamic
         * 
         * @param inv_r inverse distance between the two sites
         * @param a frozen exponent of the first site
         * @param b frozen exponent of the second site
         * @param exp_ar pre-computed exp(-a*r) using the first exponent
         * @param exp_br pre-computed exp(-a*r) using the secpnd exponent
         * @param dEdR The derivative of the energy w.r.t. the distance
         * @return The elec-elec energy in a.u.
         */
        double elec_elec_energy(const double inv_r, const double a, const double b, const double exp_ar, const double exp_br, double &dEdR);
        /**
         * @brief Compute the energy between a nuclei and an electron density
         * 
         * @param inv_r inverse distance between the two sites
         * @param a frozen exponent of the density site
         * @param exp_ar pre-computed exp(-a*r) using the density exponent
         * @param dEdR The derivative of the energy w.r.t. the distance
         * @return The nuclei-density energy in q.u. 
         */
        double elec_nuclei_energy(const double inv_r, const double a, const double exp_ar, double &dEdR);

        


        //****************************************************************
        //                      Miscellaneous
        //****************************************************************
        /**
         * @brief Get the delta-rho Coumomb J-matrix
         * 
         * @return 1D array of the symmetric matrix elements. Contains NxN elements
         */
        vec_d get_rho_coulomb_mat();
        /**
         * @brief Get the electroc potentia vector (linear term) 
         * 
         * @return vec_d Nx1 array
         */
        vec_d get_rho_pot_vec();
        /**
         * @brief Get the delta rhos of each site that are computed after minimizing the energy
         * 
         * @return vec_d Nx1 array
         */
        vec_d get_delta_rho();
        /**
         * @brief Not implimented yet! Return the total wall time used for computations
         * 
         * @param calculate_forces 
         */
        double get_total_time();
        /**
         * @brief Set wether or not forces should also be computed. If set to false, some
         * computation time can be saved
         * 
         * @param calculate_forces If true, then forces will be computed
         */
        void set_calc_forces(bool calculate_forces);
        /**
         * @brief Print site parameters to the console.
         * 
         * @param message string: A comment message to be printed
         * @param param_name string: the name of the parameter to be printed
         */
        void print_params(const std::string message, const std::string param_name);
        /**
         * @brief Get the name of the parameters that can be return
         * 
         * @return std::vector<std::string> A list of possible parameters to enquire about
         */
        std::vector<std::string> get_param_names();
        /**
         * @brief Get an array of parameters by it's name
         * 
         * @param param_name string: the name of the parameter to return
         * @return vec_d The value of the parameter for each site
         */
        vec_d get_params_by_name(const std::string param_name);
        /**
         * @brief Perform the dot product between two coordinates
         * 
         * @param coords A 3N array of coordinates for ALL sites
         * @param i index of the first first site
         * @param j index of the second site
         * @return the dot product
         */
        static double dot3Vec(const vec_d &coords, int i, int j);
        /**
         * @brief The quadratic A-matrix that was used in the last energy minimization. This includes
         * The constraint matrix in it's off-diagonal elements
         */
        vec_d A_mat_save;
        /**
         * @brief The linear B-vector used in the last energy minimization.
         * 
         */
        vec_d B_vec_save;

        enum Solver{
            Global = 0,
            DivideConquer = 1
        };

        void set_solver(Solver solver_in);


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
        int dampening_type;
        

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
        bool has_ext_field;
        Vec3 ext_field;
        vec_d ext_field_potential;

        std::vector<int> site_frag_ids;
        std::vector<std::vector<int> > site_frag_ids_partitioned;

        Timer timers;
        Solver solver;
        DivideAndConquer divide_and_conquer;
};



#endif  // FLUC_DENS_H