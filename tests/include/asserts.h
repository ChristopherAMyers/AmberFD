
#include <stdexcept>
#include <stdio.h>
#include <sstream>

#ifndef TEST_ASSERTIONS
#define TEST_ASSERTIONS


void assert_equal_tol(double found, double expected, double tol, bool print_res=true)
{
    double denom = std::min(abs(found), abs(expected));
    double diff = std::abs(found - expected);
    double pct_diff = diff/denom;
    if (print_res)
    {
        printf("\n    ASSERTION TEST    \n");
        printf(" found     : %.16f\n", found);
        printf(" expected  : %.16f\n", expected);
        printf(" tolerance : %.16f\n", tol);
        printf(" pct diff  : %.16f\n", pct_diff);
    }
    if (std::abs(pct_diff) > tol)
    {
        std::stringstream message;
        message << "ASSERTION FAILED; EXPECTED " << expected << ", FOUND " << found;
        //printf(" %s\n", message.str().c_str());
        throw std::runtime_error(message.str());
    }
}

#endif // TEST_ASSERTIONS