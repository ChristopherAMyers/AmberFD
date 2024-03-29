set(TEST_SOURCE ${PROJECT_SOURCE_DIR}/tests)

include_directories(${TEST_SOURCE}/include)
message(" TESTS WILL USE: Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
set(TEST_BUILD ${TEST_SOURCE}/build)


add_executable(test_UU_energy test_UU_energy.cpp)
add_executable(test_fluc_terms test_fluc_terms.cpp)
set_target_properties(test_UU_energy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BUILD})
set_target_properties(test_fluc_terms PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BUILD})
target_link_libraries(test_UU_energy PRIVATE AmberFD)
target_link_libraries(test_fluc_terms PRIVATE AmberFD)

add_test(NAME UU_energy  WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${TEST_BUILD}/test_UU_energy)
add_test(NAME fluc_terms  WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${TEST_BUILD}/test_fluc_terms)
add_test(NAME swig_import  WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${Python3_EXECUTABLE} ${TEST_SOURCE}/test_swig_import.py)
add_test(NAME UU_energy_py WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${Python3_EXECUTABLE} ${TEST_SOURCE}/test_UU_energy.py)
add_test(NAME DispPauli_energy_py WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${Python3_EXECUTABLE} ${TEST_SOURCE}/test_disp_energy.py)
add_test(NAME AmberFD_energy_py WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${Python3_EXECUTABLE} ${TEST_SOURCE}/test_AmberFD.py)
add_test(NAME forces_py WORKING_DIRECTORY ${TEST_SOURCE} COMMAND ${Python3_EXECUTABLE} ${TEST_SOURCE}/test_forces.py)