%{
    #define SWIG_FILE_WITH_INIT
    #include "include/DispersionPauli.h"
%}

%apply (int* IN_ARRAY1, int DIM1)    {(int *nuclei, int len1)}
%apply (double* IN_ARRAY1, int DIM1) {(double *exponents, int len2),
                                      (double *coeff, int len3),
                                      (double* coeff_list, int len),
                                      (double* exp_list, int len)}

%init %{
    import_array();
%}
%include "include/DispersionPauli.h"
%include "include/common.h"

%exception DispersionPauli::DispersionPauli {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%extend DispersionPauli
{
    DispersionPauli(int *nuclei, int len1, double *exponents, int len2, double *coeff, int len3)
    {
        if ((len1 != len2) || (len1 != len3))
            PyErr_Format(PyExc_ValueError,
                        "Arrays of lengths (%d,%d,%d) given",
                        len1, len2, len3);
        DispersionPauli* dispPauli = new DispersionPauli(len1, nuclei, exponents, coeff);
        return dispPauli;
    }
}