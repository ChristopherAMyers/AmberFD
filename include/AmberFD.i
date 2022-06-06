%module _AmberFD

//%feature("autodoc", "1");
%feature("doxygen:ignore:transferfull");

%{
    #define SWIG_FILE_WITH_INIT
    #include "include/AmberFD.h"
    //#include "include/Vec3.h"
%}

%include "numpy.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_map.i"
%include "std_set.i"
%include "std_string.i"
%include <std_shared_ptr.i>
%include <typemaps.i>




//  TODO: Test for performance hits
%include exception.i       
%exception {
    try {
        $action
    }
    catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError,std::string(e.what()).c_str());
        //PyErr_SetString(PyExc_IndexError,"index out-of-bounds");
        return NULL;
    }
}



namespace std {
    %template(MapID) map<int,double>;
    %template(PairII) pair<int,int>;
    %template(VectorD) vector<double>;
    %template(VectorDD) vector< vector<double> >;
    %template(VectorI) vector<int>;
    %template(VectorII) vector < vector<int> >;
    %template(VectorPairII) vector< pair<int,int> >;
    %template(Matrix) vector<Vec3>;
    %template(SetI) set<int>;
};

%init %{
    import_array();
%}


%shared_ptr(FlucDens);
%shared_ptr(DispersionPauli);
%include "include/FlucDens.i"
%include "include/DispersionPauli.i"
%include "include/AmberFD.h"
%include "include/common.h"
%include "include/Energies.h"
