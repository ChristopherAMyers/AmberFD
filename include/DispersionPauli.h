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

/**
 * This class impliments two forces, a repulsive Pauli force and an attractive
 * dispersion force. The sum of these two constitute a force similar to a Lennard-Jones
 * force used in many MD simulations or a van der Waals force used in the AMOEBA model.
 * To use it, simple create a DispersionPauli object by specifying the number of 
 * sites and the parameters for all sites. These parameters can be altered
 * after the creation of the force object.
 
 * 
 * Next, set the global parameters for dispersion dampening, which employs the same 
 * functional form used in Grimme et al https://doi.org/10.1002/jcc.21759 for the 
 * r^6 attractive term. If no dampening is desired, set both a1 and a2 parameters to zero. 
 * As a result, vdw Radii will not affect the dispersion energy, and only C6 coefficeients 
 * will matter. Both C6 Nd vds Radii are determined based on the nuclei of each site, and
 * all sites with the same nuclei will use the same parameters.
 * 
 * For the interaction betweena pair of sites, a geometric mixing is employed for the C6
 * coefficients while the two vdw radii are simply added together. For the Pauli
 * repulsion, exponents are arithmetically averaged while the radii are converted to 
 * coefficient pre-factors that are then multipled by each other.
 * 
 * @brief To construct a new force object, define the object with
 * 
 * @param num_sites The number of sites to create
 * @param nuclei The integer nuclei of each site
 * @param exponents The exponents used for the Pauli repulsion
 * @param radii The Pauli radii used for repulsion
 * 
 */
class DispersionPauli {

    public:
        /**
         * @brief Construct a new DispersionPauli force object
         * 
         * @param num_sites The number of sites to create
         * @param nuclei The integer nuclei of each site
         * @param exponents The exponents used for the Pauli repulsion
         * @param radii The Pauli radii used for repulsion
         */
        DispersionPauli(const int num_sites,
                        const int* nuclei,
                        const double* exponents,
                        const double* radii);
        ~DispersionPauli();

        //****************************************************************
        //                       Model Parameters
        //****************************************************************
        /**
         * @brief Set the dispersion parameters assording to 
         * Grimme et alhttps://doi.org/10.1002/jcc.21759
         * 
         * @param s6 Coefficient that scales the entire energy
         * @param a1 multiple of the vdw radii R_ab
         * @param a2 additive to a1*R_ab
         */
        void set_dispersion_params(double s6, double a1, double a2);
        /**
         * @brief Get the dispersion parameters
         * 
         * @param s6 Coefficient that scales the entire energy
         * @param a1 multiple of the vdw radii R_ab
         * @param a2 additive to a1*R_ab
         */
        void get_dispersion_params(double &s6, double &a1, double &a2);
        /**
         * @brief Set the vdw radii by the nuclei type
         * 
         * @param nuclei2radiiMap mapping from nuclei integers to vdw Radii (in bohr)
         */
        void set_vdw_radii(map_id nuclei2radiiMap);
        /**
         * @brief Get the vdw radii mapping by the nuclei type
         * 
         * @return mapping from nuclei integers to vdw Radii (in bohr) 
         */
        map_id get_vdw_radii_map();
        /**
         * @brief Get the vdw radii used for all sites
         * 
         * @return array of vdw radii in bohr 
         */
        vec_d get_vdw_radii();
        /**
         * @brief Set the C6 coefficients by the nuclei type
         * 
         * @param nucleiToC6Map mapping from nuclei integers to C6 coeff
         */
        void set_C6_map(map_id nucleiToC6Map);
        /**
         * @brief Get the C6 coefficient mapping by nuclei type
         * 
         * @return mapping fron nuclei integers to C6 coeff
         */
        map_id get_C6_map();
        /**
         * @brief Get the C6 coefficients used for all sites 
         * 
         * @return array of C6 coefficients
         */
        vec_d get_C6_coeff();
        /**
         * @brief Set the Pauli radii for all sites
         * 
         * @param radii_list list of Pauli radii
         */
        void set_pauli_radii(vec_d radii_list);
        /**
         * @brief Set the Pauli radii for an individual site
         * 
         * @param index the index of the site to change
         * @param radii the Pauli radii to assign
         */
        void set_pauli_radii(int index, double radii);
        /**
         * @brief Get the Pauli radii for all sites
         * 
         * @return array of Pauli radii in bohr
         */
        vec_d get_pauli_radii();
        /**
         * @brief Set the Pauli exponents for all sites
         * 
         * @param exp_list array of exponents to use
         */
        void set_pauli_exp(vec_d exp_list);
        /**
         * @brief Set the Pauli exponent for a particular site
         * 
         * @param index The index of the site to change the exponent of
         * @param exponent The value of the exponent to use on this site
         */
        void set_pauli_exp(int index, double exponent);
        /**
         * @brief Get the Pauli exponents used on all sites in the system
         * 
         * @return vec_d Array of Pauli exponents (in inverse bohr)
         */
        vec_d get_pauli_exp();
        /**
         * @brief Get the number sites
         * 
         * @return int 
         */
        int get_num_sites();

        //****************************************************************
        //                      Exclusions
        //****************************************************************
        /**
         * @brief Create exclusions based on wht atoms are bonded to each other
         * 
         * @param bonds A list of index pairs in which each pair declares a bond between
         * those two sites
         * @param bond_cutoff Any site this many bonds away are excluded in the interaction
         * calculation 
         */
        void create_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, int bond_cutoff);
        /**
         * @brief Create exclusions from fragments. A fragment is a collection of indicies in which
         * all sites will not interact with each other.
         * 
         * @param frag_idx An array of indicies
         */
        void create_exclusions_from_fragment(const vec_i frag_idx);
        /**
         * @brief Add an exclusion between a pair of sites
         * 
         * @param i The index of the first site
         * @param j The index of the second site
         */
        void add_exclusion(const int i, const int j);
        /**
         * @brief Get a list of sites that are excluded from interacting with the i-th site
         * 
         * @param particle1 The site to inquire about exclusions
         * @return std::set<int> a list of site indicies that do not interact with particle1
         */
        std::set<int> get_exclusions(const int particle1) const;

        //****************************************************************
        //                  Periodic Boundary Condtions
        //****************************************************************
        /**
         * @brief Turn on or off periodic boundary conditions
         * 
         * @param is_periodic If true, use PBC, if false, do not use them
         */
        void set_use_PBC(bool is_periodic);
        /**
         * @brief Set the periodic boundaries
         * 
         * @param is_periodic If true, use PBC, if false, do not use them
         * @param x The x-direction of the periodic box
         * @param y The y-direction of the periodic box
         * @param z The z-direction of the periodic box
         */
        void set_use_PBC(bool is_periodic, const double x, const double y, const double z);
        /**
         * @brief Inquire if PBC is used or not
         * 
         * @return true PBC is being used
         * @return false PBC is not being used
         */
        bool get_use_PBC();
        
        //****************************************************************
        //                          Energies
        //****************************************************************
        /**
         * @brief Clear energy components before a calculation
         * 
         */
        void initialize();
        /**
         * @brief Calculate the total interaction energy of all sites
         * 
         * @param coords A 1D array of positions of the sites
         * @return The interaction energy
         */
        double calc_energy(const vec_d &coords);
        /**
         * @brief Compute the interaction energy of a single pair of sites
         * 
         * @param pos  A 1D array of positions of the sites
         * @param deltaR A distance object for the pair of sites
         * @param i Index of the first site
         * @param j Index of the second site
         * @param energies The Energies object that will be updated
         * @param forces The forces array that will be updated
         * @return double 
         */
        double calc_one_pair(const vec_d &pos, DeltaR &deltaR, int i, int j, Energies& energies, std::vector<Vec3> &forces);
        /**
         * @brief Compute the interaction energy of a single pair of sites
         * 
         * @param pos A 1D array of positions of the sites
         * @param i Index of the first site
         * @param j Index of the second site
         * @param forces The forces array that will be updated
         * @return Energies that were computed for the pair
         */
        Energies calc_one_pair(const vec_d &pos, int i, int j, std::vector<Vec3> &forces);
        /**
         * @brief Return the total Pauli repulsion energy for the system
         * 
         * @return double 
         */
        double get_pauli_energy();
        /**
         * @brief Return the total dispersion energy for the system
         * 
         * @return double 
         */
        double get_disp_energy();
        /**
         * @brief Return the forces exerted on each site after an energy computation
         * 
         * @return std::vector<vec_d> An Nx3 array of forces
         */
        std::vector<vec_d> get_forces();

        //****************************************************************
        //                      Miscellaneous
        //****************************************************************
        /**
         * @brief Set the use of secondary radii for Pauli repulsion. This is a backup
         * to ensure that two sites never get too close to each other without an enormous
         * energetic cost. The distances used are chosed to be a close to a third of the 
         * vdw radii so that the energy cost will be minimal at normal distances. This 
         * can redude performance somewhat.
         * 
         * @param use_radii If true, the backup radii will alos be used 
         */
        void set_use_secondary_radii(bool use_radii=true);
        /**
         * @brief Set the use two site repulsion (experimental only)
         * 
         * @param on_off Set to true to enable
         */
        void set_use_two_site_repulsion(bool on_off);
        /**
         * @brief Set the two site repulsion distance
         * 
         * @param vertical_dist The virtical distance out of the place
         */
        void set_two_site_distance(double vertical_dist);
        /**
         * @brief Create the two-site repulsion sites
         * 
         * @param vertical_dist The virtical distance out of the place
         * @param bonds A list of all bonds
         */
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

        std::vector<Vec3> self_forces;
        //std::vector<vec_d> self_forces;
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