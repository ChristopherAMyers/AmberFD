#include "DispersionPauli.h"

DispersionPauli::DispersionPauli(const int num_sites, const int* nuclei_in, const double* exponents_in, const double* coeff_in)
{
    n_sites = num_sites;

    //  copy over parameters
    nuclei.assign(nuclei_in, nuclei_in + n_sites);
    exponents.assign(exponents_in, exponents_in + n_sites);
    coeff.assign(coeff_in, coeff_in + n_sites);

    //  initialize exclusions
    exclusions.resize(num_sites);
}

DispersionPauli::~DispersionPauli(){}