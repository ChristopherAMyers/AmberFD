
project(FlucCharge)
cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# include and sources
set(INCL ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(include)
include_directories(/usr/include)
include_directories(/usr/include/x86_64-linux-gnu)
file(GLOB_RECURSE SOURCES src/*.cpp)

# OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    add_compile_options(-O3 -fopenmp)
    message("Setting OpenMP Flags")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# python libraries
find_package(Python3 COMPONENTS Interpreter Development NumPy)
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

#   lapack and blas
find_package(LAPACK REQUIRED)
include_directories(${BLAS_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR})

# SWIG
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS -doxygen)
set_source_files_properties(${INCL}/AmberFD.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(${INCL}/FlucDens.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(${INCL}/swigDispersionPauli.i PROPERTIES CPLUSPLUS ON)
swig_add_library(swigAmberFD TYPE SHARED LANGUAGE python SOURCES ${INCL}/AmberFD.i)
swig_link_libraries(swigAmberFD AmberFD ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

#   print out python versions and locations
message("PYTHONLIBS_VERSION_STRING: ${PYTHONLIBS_VERSION_STRING}")
message("Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message("Python3_LIBRARIES: ${Python3_LIBRARIES}")
message("BLAS_INCLUDE_DIR: ${BLAS_INCLUDE_DIR}")
message("BLAS_LIBRARIES ${BLAS_LIBRARIES}")
message("LAPACK_INCLUDE_DIR: ${LAPACK_INCLUDE_DIR}")
message("LAPACK_LIBRARIES ${LAPACK_LIBRARIES}")
message("CMAKE_SWIG_FLAGS: ${CMAKE_SWIG_FLAGS}")
message("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
message("INCL: ${INCL}")
foreach(file ${files})
    message("test file: ${file}")
endforeach()

# compile options
add_compile_options(-fPIC -O3 -g)

# library compilation
add_library(AmberFD STATIC ${SOURCES})
target_link_libraries(AmberFD ${BLAS_LIBRARIES})
target_link_libraries(AmberFD ${LAPACK_LIBRARIES})
target_link_libraries(AmberFD lapacke)

#   tests
target_include_directories(AmberFD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
enable_testing()
add_subdirectory(tests)