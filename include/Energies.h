#ifndef ENERGIES_H_
#define ENERGIES_H_

/**
 * @brief A simple struct that keeps track of the energy components
 * throughout an computation. AmberFD, DispersionPauli, and FlucDens
 * forces may return these objects after a calculation.
 */
class Energies{
    public:
    /**
     * @brief Construct a new Energies object with all energies zeroed out.
     * 
     */
        Energies()
        {
            zero();
        }
        /**
         * @brief Construct a new Energies object (copy constructor)
         * 
         * @param Energies 
         */
        Energies(const Energies &eng)
        {
            pauli      = eng.pauli; /** Pauli repulsion */
            disp       = eng.disp;
            frz        = eng.frz;
            pol        = eng.pol;
            vct        = eng.vct;
            elec_elec  = eng.elec_elec;
            elec_nuc   = eng.elec_nuc;
            nuc_nuc    = eng.nuc_nuc;
            pauli_wall = eng.pauli_wall;
            frz_ext    = eng.frz_ext;
        }
        /**
         * @brief Return the total energy. This is the sum of pauli, disp, frz, pol, vct, frz_ext, and pauli_wall
         * 
         * @return double 
         */
        double total()
        {
            return pauli + disp + frz + pol + vct + pauli_wall + frz_ext;
        }
        /**
         * @brief Zero out all the energies.
         * 
         */
        void zero()
        {
            pauli = disp = frz = pol = vct = pauli_wall = 0.0;
            elec_elec = elec_nuc = nuc_nuc = 0.0;
            frz_ext = 0.0;
        }
        /**
         * @brief Sum another Energies object into this one
         * 
         * @param eng The Energies object to be summed from
         */
        void add(const Energies &eng)
        {
            pauli      += eng.pauli;
            disp       += eng.disp;
            frz        += eng.frz;
            pol        += eng.pol;
            vct        += eng.vct;
            elec_elec  += eng.elec_elec;
            elec_nuc   += eng.elec_nuc;
            nuc_nuc    += eng.nuc_nuc;
            pauli_wall += eng.pauli_wall;
            frz_ext    += eng.frz_ext;
        }

        /**
         * @brief Addition operator
         */
        Energies& operator+=(const Energies &eng)
        {
            add(eng);
            return *this;
        }

        // Energies operator-();
        // Energies& operator-=(const Energies &rhs);
        // Energies operator*();
        // Energies& operator*=(double rhs);
        // Energies operator/();
        // Energies& operator/=(double rhs);

        /** Exponential Pauli repulsion */
        double pauli; 
        /** Dispersion energy */
        double disp;
        /** Frozen electrostatics, this is NOT the frozen energy from an energy decomposition analysis */
        double frz; 
        /** Polarization energy */
        double pol; 
        /** Charge transfer energy */
        double vct; 
        /** Electron-Electron Coulomb repulsion*/
        double elec_elec; 
        /** Electron-Nuclei attraction */
        double elec_nuc; 
        /** Nuclear-Nuclear repulsion */
        double nuc_nuc; 
        /** Secondary Pauli "wall" energy */
        double pauli_wall; 
        /** Frozen electrostatics with external field. Dynamic electron density interacting with the field is already accounted for in the polarization energy*/
        double frz_ext; 

    private:   
        void multiply(double value)
        {
            pauli      *= value;
            disp       *= value;
            frz        *= value;
            pol        *= value;
            vct        *= value;
            elec_elec  *= value;
            elec_nuc   *= value;
            nuc_nuc    *= value;
            pauli_wall *= value;
            frz_ext    *= value;
        }
        void subtract(const Energies &eng)
        {
            pauli      -= eng.pauli;
            disp       -= eng.disp;
            frz        -= eng.frz;
            pol        -= eng.pol;
            vct        -= eng.vct;
            elec_elec  -= eng.elec_elec;
            elec_nuc   -= eng.elec_nuc;
            nuc_nuc    -= eng.nuc_nuc;
            pauli_wall -= eng.pauli_wall;
            frz_ext    -= eng.frz_ext;
        }
};

#endif // ENERGIES_H_