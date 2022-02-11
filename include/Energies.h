#ifndef ENERGIES_H_
#define ENERGIES_H_

class Energies{
    public:
        Energies()
        {
            zero();
        }
        Energies(const Energies &eng)
        {
            pauli      = eng.pauli;
            disp       = eng.disp;
            frz        = eng.frz;
            pol        = eng.pol;
            vct        = eng.vct;
            elec_elec  = eng.elec_elec;
            elec_nuc   = eng.elec_nuc;
            nuc_nuc    = eng.nuc_nuc;
            pauli_wall = eng.pauli_wall;
        }

        double total()
        {
            return pauli + disp + frz + pol + vct + pauli_wall;
        }
        void zero()
        {
            pauli = disp = frz = pol = vct = pauli_wall = 0.0;
            elec_elec = elec_nuc = nuc_nuc = 0.0;
        }
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
        }

        // Energies operator+();
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

        double pauli, disp, frz, pol, vct;
        double elec_elec, elec_nuc, nuc_nuc;
        double pauli_wall;

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
        }
};

#endif // ENERGIES_H_