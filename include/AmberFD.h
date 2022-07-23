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
#include <fstream>
#include <cstdio>
#include <omp.h>

#include "common.h"
#include "FlucDens.h"
#include "DispersionPauli.h"
#include "fileReader.h"

#ifndef AMBERFD_H
#define AMBERFD_H

/**  
 * Small structure used to assign per-site force field parameters.
 * This is passed into AmberFD when adding a new particle.
 */
class ParticleInfo{
    public:
        /**
         * @brief Construct a new Particle Info object
         * 
         * @param nuclei The integer nuclei of the particle. A value of zero indicates it's a virtual particle
         */
        ParticleInfo(int nuclei)
        {
            this->nuclei = nuclei;
            frz_chg = frz_exp = dyn_exp = pauli_exp =  pauli_radii = 0.0;
        }
        /** Integer value for the nuclei of the particle. A value of zero indicates it's a virtual particle */
        int nuclei;
        /** Frozen charge of the particle, can be positive or negative */
        double frz_chg;
        /** Frozen Slater exponent of the particle, must be a positive value */
        double frz_exp;
        /** Dynamic Slater exponent of the particle, must be a positive value */
        double dyn_exp;
        /** Pauli exponent of the particle, must be a positive value */
        double pauli_exp;
        /** Pauli radii of the particle, must be a positive value */
        double pauli_radii;

        // ParticleInfo(int nuclei, double frz_chg, double frz_exp, double dyn_exp, double pauli_exp, double pauli_radii):
        // nuclei(nuclei), frz_chg(frz_chg), frz_exp(frz_exp), dyn_exp(dyn_exp), pauli_exp(pauli_exp), pauli_radii(pauli_radii)
        // {}
};

/**
 * AmberFD is an umbrella function thatcombines FlucDens force 
 * and DispersionPauli force settings into one accessable class. This
 * class does not calculate forces itslf, but is instead used create and access the 
 * fluctuating density and dispersion-pauli system and set it's options. 
 * 
 * To use this class,  first create an AmberFD object with the number of particles to be added to 
 * the system. The, add each particle to the system by calling add_particle(),
 * passing in a ParticleInfo object. Once all particles have been added, create
 * the fluctuating density force and dispersion-pauli force by calling create_fluc_dens_force()
 * and create_disp_pauli_force(), respectively. Each of these forces will then be generated,
 * and a pointer to each force can be obtained by calling get_fluc_dens_force() or 
 * get_disp_pauli_force().
 * 
 * Once the forces have been generated, adding new particles will have no effect on the system.
 * Calling create_XXX_force() will simply return a pointer to the already created forces.
 * If periodic boundary conditions are desired, then call set_use_PBC() after all particles have
 * been added.
 *  
 */

class AmberFD{

    public:
        /**
         * Construct a new AmberFD force
         */
        AmberFD();
        /**
         * Construct a new AmberFD force with the number of particles already known.
         * This can help speed up the system creation for very large numbers of particles.
         * More particles can be added than n_sites, but the memory will not be preallocated.
         * 
         * @param n_sites the estimated number of particles to create a system for 
         */
        AmberFD(const int n_sites);
        ~AmberFD();
        /**
         * Add a new particle to the system with the properties specified in it's ParticleInfo
         * 
         * @param index 
         * @param parameters 
         * @return index of the particle added
         */
        int add_particle(int index, ParticleInfo parameters);
        /**
         * Set the indicies that define a fragment within the entire molecular system.
         * 
         * @param frag_idx 
         */
        void add_fragment(const vec_i frag_idx);
        /**
         * Get the number of particles in the system
         * 
         * @return int 
         */
        int get_num_particles()
        {   return n_sites; }
        /**
         * Return the forces exerted on each particle. This will only have an effect
         * if calc_energy_forces has been called forst with the updated positions. Else, this
         * will return the forces most recently computed.
         * 
         * @return std::vector<vec_d> 
         */
        std::vector<vec_d> get_forces();
        /**
         * Compute the energy and forces on the AmberFD system at the current positions.
         * 
         * @param positions 1-D array of positions (in a.u.)
         * @return Energies object
         */
        Energies calc_energy_forces(const vec_d &positions);
        /**
         * Get a pointer to the FlucDens force object being used
         * 
         * @param create_if_null create the force if it hasn't been created yet.
         * @return std::shared_ptr<FlucDens> 
         */
        std::shared_ptr<FlucDens> get_fluc_dens_force(bool create_if_null=false);
        /**
         * @brief Get a pointer to the DispersionPauli force object being used
         * 
         * @param create_if_null create the force if it hasn't been created yet.
         * @return std::shared_ptr<DispersionPauli> 
         */
        std::shared_ptr<DispersionPauli> get_disp_pauli_force(bool create_if_null=false);
        /**
         * @brief Create a FlucDens force object. Adding particles to the system
         *        After creating this force will have no effect.
         * 
         * @return std::shared_ptr<FlucDens> a pointer to the created force
         */
        std::shared_ptr<FlucDens> create_fluc_dens_force();
        /**
         * @brief Create a DispersionPauli force object. Adding particles to the system
         *        After creating this force will have no effect.
         * 
         * @return std::shared_ptr<DispersionPauli> a pointer to the created force
         */
        std::shared_ptr<DispersionPauli> create_disp_pauli_force();
        /**
         * @brief Get the a mapping from the OpenMM particle index to
         *        the AmberFD particle index.
         * 
         * @return std::map<int, int> 
         */
        std::map<int, int> get_index_mapping();
        /**
         * @brief Calculate the interaction between a single pair of particles. Since
         *        Polarization is a many-body force, it will not be calculated here.
         * 
         * @param positions The full array of positions of the system that would normally be passed to calc_energy_forces()
         * @param i The index of the first particle
         * @param j The index of the second particle
         * @return Energies The energies of the pair
         */
        Energies calc_one_pair(const vec_d &positions, int i, int j);
        /**
         * @brief Get the distance object between two pairs of particles.
         * 
         * @param positions The full array of positions of the system that would normally be passed to calc_energy_forces()
         * @param i The index of the first particle
         * @param j The index of the second particle
         * @return DeltaR 
         */
        DeltaR getDeltaR(const vec_d &positions, int i, int j);
        /**
         * @brief Turn on or off the use of periodic boundary conditions.
         * 
         * @param is_periodic 
         */
        void set_use_PBC(const bool is_periodic);
        /**
         * @brief Turn on or off the use of periodic boundary conditions.
         * 
         * @param is_periodic 
         * @param x The x-direction box length
         * @param y The y-direction box length
         * @param z The z-direction box length
         */
        void set_use_PBC(const bool is_periodic, const double x, const double y, const double z);
        /**
         * @brief Determine of the periodic doundary contitions are used or not
         * 
         * @return true 
         * @return false 
         */
        bool get_use_PBC();
        /**
         * @brief Serialize force FlucDens and DisperionPauli forces to a file
         * 
         * @param file_loc String of the file location to write data to.
         */
        void dump_to_file(std::string file_loc);
        /**
         * @brief Load in serialized data previously created with dump_to_file()
         * 
         * @param file_loc String of the file location to load.
         */
        void load_from_file(std::string file_loc);
        /**
         * @brief Set the number of threads to be used in force computations vis OpenMP.
         *        If set less than or equal to one, then a serial calculation will take place.
         * 
         * @param n_threads Number of threads to use
         */
        void set_threads(int n_threads);
        /**
         * @brief Get the cumulative wall time used to compute forces and energies so far.
         * 
         * @return Wall time in second
         */
        double get_parallel_time();

        ///////////////   NEW
        void print_timings();

    private:
        vec_i nuclei;
        vec_d frz_chg;
        vec_d frz_exp;
        vec_d dyn_exp;
        vec_d pauli_exp;
        vec_d pauli_radii;
        std::map<int, int> omm_to_index;
        std::map<int, int> index_to_omm;
        int n_sites;
        void initialize();

        std::shared_ptr<FlucDens> flucDens;
        std::shared_ptr<DispersionPauli> dispersionPauli;

        //void zero_energies();
        Energies total_energies;
        double E_total;
        std::vector<Vec3> self_forces;
        Periodicity periodicity;

        //  multithreading
        Energies calc_threaded_energy(const vec_d &positions);
        std::vector<std::vector<Vec3>> thread_forces;
        std::vector<Energies> thread_energies;
        double parallel_time;

        Timer timers;
};



#endif // AMBERFD_H