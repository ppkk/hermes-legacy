project(hartree-fock-neweigen)
add_executable(${PROJECT_NAME} main.cpp)
include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
