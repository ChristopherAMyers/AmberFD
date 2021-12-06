//%module FlucDens
%{
    #define SWIG_FILE_WITH_INIT
    #include "include/FlucDens.h"
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* arr1, int len1),
                                      (double* arr2, int len2)}
%apply (double* IN_ARRAY1, int DIM1) {(double *frozen_chg, int len1),
                                      (double *nuclei, int len2),
                                      (double *frozen_exp, int len3),
                                      (double *dynamic_exp, int len4)}

%include "include/FlucDens.h"
%include "include/common.h"

%exception FlucDens::FlucDens {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%extend FlucDens
{
    FlucDens(double *frozen_chg, int len1, double *nuclei, int len2, double *frozen_exp, int len3, double *dynamic_exp, int len4)
    {
        if ((len1 != len2) || (len1 != len3) || (len1 != len4))
            PyErr_Format(PyExc_ValueError,
                        "Arrays of lengths (%d,%d,%d,%d) given",
                        len1, len2, len3, len4);
        FlucDens* fluc = new FlucDens(len1, frozen_chg, nuclei, frozen_exp, dynamic_exp);
        return fluc;
    }
}