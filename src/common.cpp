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