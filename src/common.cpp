#include "common.h"


void Nonbonded::calc_dR(const vec_d &coords, int i, int j, double* deltaR)
{
    deltaR[XIdx] = coords[i]     - coords[j];
    deltaR[YIdx] = coords[i + 1] - coords[j + 1];
    deltaR[ZIdx] = coords[i + 2] - coords[j + 2];
    deltaR[R2Idx] = dot3(deltaR, deltaR);
    deltaR[RIdx] = sqrt(deltaR[R2Idx]);
    deltaR[RInvIdx] = 1/deltaR[RIdx];
}

double Nonbonded::dot3(const double* u, const double* v)
{
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

std::vector<std::set<int> > Nonbonded::calc_exclusions_from_bonds(const std::vector<std::pair<int, int> > bonds, const int bond_cutoff, const int n_sites)
{
    /*  determine exclusions for frozen-frozen interactions from bonds up to 
        an arbitrary bond neighbor length. Edited from OpenMM */
    std::vector<std::set<int> > exclusions;
    if (bond_cutoff < 1)
        return exclusions;
    for (auto& bond : bonds)
    {
        if (bond.first < 0 || bond.second < 0 || bond.first >= n_sites || bond.second >= n_sites)
            throw "createExclusionsFromBonds: Illegal particle index in list of bonds";
    }
    exclusions.resize(n_sites);
    std::vector<std::set<int> > bonded12(exclusions.size());
    for (auto& bond : bonds) {
        int p1 = bond.first;
        int p2 = bond.second;
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
        bonded12[p1].insert(p2);
        bonded12[p2].insert(p1);
    }
    for (int level = 0; level < bond_cutoff-1; level++) {
        std::vector<std::set<int> > currentExclusions = exclusions;
        for (int i = 0; i < (int) n_sites; i++)
            for (int j : currentExclusions[i])
                exclusions[j].insert(bonded12[i].begin(), bonded12[i].end());
    }
    return exclusions;
}