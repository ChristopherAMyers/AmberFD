#include "common.h"

void Nonbonded::add_Vec3_to_vector(std::vector<double> &vec, const Vec3 &vec3)
{
    vec[0] += vec3[0];
    vec[1] += vec3[1];
    vec[2] += vec3[2];
}

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

DeltaR::DeltaR(double *deltaR)
{
    dR = Vec3(deltaR[Nonbonded::XIdx], deltaR[Nonbonded::YIdx], deltaR[Nonbonded::ZIdx]);
    r2 = deltaR[Nonbonded::R2Idx];
    r = deltaR[Nonbonded::RIdx];
    r_inv = deltaR[Nonbonded::RInvIdx];
}
DeltaR::DeltaR(const vec_d &coords, int i, int j)
{
    getDeltaR(coords, i, j);
}
DeltaR::DeltaR(const vec_d &coords, int i, int j, Periodicity p)
{
    getDeltaR(coords, i, j, p);
}

void DeltaR::getDeltaR(const vec_d &coords, int i, int j)
{
    double dx = coords[i + 0] - coords[j + 0];
    double dy = coords[i + 1] - coords[j + 1];
    double dz = coords[i + 2] - coords[j + 2];
    dR = Vec3(dx, dy, dz);
    r2 = dR.dot(dR);
    r = sqrt(r2);
    r_inv = 1/r;
}

void DeltaR::getDeltaR(const vec_d &coords, int i, int j, Periodicity p)
{
    double dx = coords[i + 0] - coords[j + 0];
    double dy = coords[i + 1] - coords[j + 1];
    double dz = coords[i + 2] - coords[j + 2];
    dR = Vec3(dx, dy, dz);
    dR[2] -= p.box_size[2]*ceil(dR[2]*p.inv_box_size[2] - 0.5);
    dR[1] -= p.box_size[1]*ceil(dR[1]*p.inv_box_size[1] - 0.5);
    dR[0] -= p.box_size[0]*ceil(dR[0]*p.inv_box_size[0] - 0.5);
    r2 = dR.dot(dR);
    r = sqrt(r2);
    r_inv = 1/r;
    //printf(" IN getDeltaR: %.10f  %.10f  %.10f  %.10f \n", p.inv_box_size[2], dR[0], dR[1], dR[2]);
}

void DeltaR::getDeltaR(const Vec3 &ri, const Vec3 &rj, Periodicity p)
{
    dR = ri - rj;
    dR[2] -= p.box_size[2]*ceil(dR[2]*p.inv_box_size[2] - 0.5);
    dR[1] -= p.box_size[1]*ceil(dR[1]*p.inv_box_size[1] - 0.5);
    dR[0] -= p.box_size[0]*ceil(dR[0]*p.inv_box_size[0] - 0.5);
    r2 = dR.dot(dR);
    r = sqrt(r2);
    r_inv = 1/r;
}

void DeltaR::getDeltaR(const Vec3 &ri, const Vec3 &rj)
{
    dR = ri - rj;
    r2 = dR.dot(dR);
    r = sqrt(r2);
    r_inv = 1/r;
}

void DeltaR::get_pointer(double *deltaR)
{
    //double deltaR[Nonbonded::RMaxIdx];
    deltaR[Nonbonded::XIdx] = dR[0];
    deltaR[Nonbonded::YIdx] = dR[1];
    deltaR[Nonbonded::ZIdx] = dR[2];
    deltaR[Nonbonded::R2Idx] = r2;
    deltaR[Nonbonded::RIdx] = r;
    deltaR[Nonbonded::RInvIdx] = r_inv;
}

Energies::Energies()
{
    zero();
}

void Energies::add(Energies &eng)
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

void Energies::zero()
{
    pauli = disp = frz = pol = vct = pauli_wall = 0.0;
    elec_elec = elec_nuc = nuc_nuc = 0.0;
}

double Energies::total()
{
     return pauli + disp + frz + pol + vct + pauli_wall;
}


void Periodicity::set(bool is_periodic_in, const double x, const double y, const double z)
{
    is_periodic = is_periodic_in;
    box_size = Vec3(x, y, z);
    inv_box_size = Vec3(1/x, 1/y, 1/z);
    box_vectors[0] = Vec3(x, 0, 0);
    box_vectors[1] = Vec3(0, y, 0);
    box_vectors[2] = Vec3(0, 0, z);
    box_size[0] = x;
    box_size[1] = y;
    box_size[2] = z;
    inv_box_size[0] = 1/x;
    inv_box_size[1] = 1/y;
    inv_box_size[2] = 1/z;
}