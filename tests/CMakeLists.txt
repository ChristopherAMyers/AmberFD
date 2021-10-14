set(TEST_SOURCE ${PROJECT_SOURCE_DIR}/tests)
include_directories(${TEST_SOURCE}/include)
include_directories(${TEST_SOURCE})
set(TEST_BUILD ${TEST_SOURCE}/build)


add_executable(test_UU_energy test_UU_energy.cpp)
set_target_properties(test_UU_energy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BUILD})
target_link_libraries(test_UU_energy FlucDens)

add_test(NAME UU_energy    WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${TEST_BUILD}/test_UU_energy)
add_test(NAME UU_energy_py WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${Python3_EXECUTABLE} ${TEST_SOURCE}/test_UU_energy.py)
add_test(NAME swig_import  WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${Python3_EXECUTABLE} ${TEST_SOURCE}/test_swig_import.py)