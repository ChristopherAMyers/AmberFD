%module AmberFD

%{
    #define SWIG_FILE_WITH_INIT
    #include "include/AmberFD.h"
%}

%include "numpy.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_map.i"
%include "std_string.i"

namespace std {
    %template(MapID) map<int,double>;
    %template(PairII) pair<int,int>;
    %template(VectorD) vector<double>;
    %template(VectorDD) vector< vector<double> >;
    %template(VectorI) vector<int>;
    %template(VectorII) vector < vector<int> >;
    %template(VectorPairII) vector< pair<int,int> >;
};

%init %{
    import_array();
%}
%include "include/FlucDens.i"
%include "include/DispersionPauli.i"
%include "include/AmberFD.h"
%include "include/common.h"
