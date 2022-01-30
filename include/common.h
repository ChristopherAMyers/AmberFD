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

class Periodicity{
    public:
        Periodicity() : is_periodic(false){};
        bool is_periodic;
        Vec3 box_vectors[3];
        Vec3 box_size;
        Vec3 inv_box_size;
        void set(const bool is_periodic, const double x, const double y, const double z);
};

class DeltaR{
    public:
        DeltaR() = default;
        DeltaR(double *deltaR);
        DeltaR(const vec_d &coords, int i, int j);
        DeltaR(const vec_d &coords, int i, int j, Periodicity period);
        void getDeltaR(const vec_d &coords, int i, int j);
        void getDeltaR(const vec_d &coords, int i, int j, Periodicity period);
        Vec3 dR;
        double r, r_inv, r2;
        void get_pointer(double *deltaR);
};

class Energies{
    public:
        double pauli, disp, frz, pol, vct;
        double elec_elec, elec_nuc, nuc_nuc;
        double pauli_wall;
        Energies();
        double total();
        void zero();
        void add(Energies &eng);
};



#endif // TYPEDEFS_H