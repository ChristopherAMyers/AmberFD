#include <vector>
#include <math.h>
#include <set>
#include <map>
#include <utility>

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#define ANG2BOHR 1.8897259886

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

        static void calc_dR(const vec_d &coords, int i, int j, double* deltaR);
        static double dot3(const double*u, const double* v);

        static std::vector<std::set<int> > calc_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, const int bond_cutoff, const int n_sites);

};

#endif // TYPEDEFS_H