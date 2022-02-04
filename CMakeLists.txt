
project(FlucCharge)
cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#SET(Python3_INCLUDE_DIRS /home/cmyers/miniconda3/include/python3.7m)
#SET(Python3_EXECUTABLE /home/cmyers/miniconda3/bin/python3)
#SET(Python3_LIBRARIES /home/cmyers/miniconda3/lib/libpython3.7m.so)

find_package(Python3 COMPONENTS Interpreter Development NumPy)
find_package(SWIG REQUIRED)
#find_package(OpenMP REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS -py3)

#   lapack and blas
#find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#set(LAPACK_INCLUDE_DIR /usr/include)

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
set(INCL ${CMAKE_CURRENT_SOURCE_DIR}/include)
message("INCL: ${INCL}")

include_directories(/usr/include)
include_directories(/usr/include/x86_64-linux-gnu)
include_directories(${BLAS_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR})
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(include)
set_source_files_properties(${INCL}/AmberFD.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(${INCL}/FlucDens.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(${INCL}/swigDispersionPauli.i PROPERTIES CPLUSPLUS ON)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

file(GLOB_RECURSE SOURCES src/*.cpp)
foreach(file ${files})
    message("test file: ${file}")
endforeach()

add_compile_options(-fPIC -O2 -pg)
if (OPENMP_FOUND)
    message("Setting OpenMP Flags")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(AmberFD STATIC ${SOURCES})
target_link_libraries(AmberFD ${BLAS_LIBRARIES})
target_link_libraries(AmberFD ${LAPACK_LIBRARIES})
#target_link_libraries(AmberFD /network/rit/lab/ChenRNALab/bin/OpenBLAS-0.3.6/lib/libopenblas.so)
#target_link_libraries(AmberFD /usr/lib/x86_64-linux-gnu/liblapacke.so)
target_link_libraries(AmberFD lapacke)

#   needed for tests
target_include_directories(AmberFD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

swig_add_library(swigAmberFD TYPE SHARED LANGUAGE python SOURCES ${INCL}/AmberFD.i)
swig_link_libraries(swigAmberFD AmberFD ${PYTHON_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

enable_testing()
add_subdirectory(tests)