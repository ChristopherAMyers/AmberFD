#include <vector>
#include <math.h>
#include <set>
#include <map>
#include <utility>
#include "Vec3.h"

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#define ANG2BOHR 1.8897259886
#define AU_2_KJ_PER_MOL 2625.5009

typedef std::vector<double> vec_d;
typedef std::vector<float> vec_f;
typedef std::vector<int> vec_i;
typedef std::map<int, double> map_id;

class Nonbonded{

    public:

        //  indicies that reference dR information
        static const int XIdx = 0;
        static const int YIdx = 1;
        static const int ZIdx = 2;
        static const int RIdx = 3;
        static const int R2Idx = 4;
        static const int RInvIdx = 5;
        static const int RMaxIdx = 6;

        static void add_Vec3_to_vector(std::vector<double> &vec, const Vec3 &vec3);

        static void calc_dR(const vec_d &coords, int i, int j, double* deltaR);
        static double dot3(const double*u, const double* v);

        static std::vector<std::set<int> > calc_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, const int bond_cutoff, const int n_sites);
};

class DeltaR{
    public:
        DeltaR() = default;
        DeltaR(double *deltaR);
        DeltaR(const vec_d &coords, int i, int j);
        Vec3 dR;
        double r, r_inv, r2;
        void get_pointer(double *deltaR);
};

class Energies{
    public:
        double pauli, disp, frz, pol, vct;
        double elec_elec, elec_nuc, nuc_nuc;
        double pauli_wall;
        Energies(){   reset_all();  }
        void reset(){   pauli = disp = frz = pol = vct = pauli_wall = 0.0; }
        double total() { return pauli + disp + frz + pol + vct + pauli_wall;   }

        void reset_all()
        { reset(); elec_elec = elec_nuc = nuc_nuc = 0.0;    }

        void add(Energies &eng)
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
};

#endif // TYPEDEFS_H