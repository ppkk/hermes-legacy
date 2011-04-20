project(magpot)

add_executable(${PROJECT_NAME} main.cpp)
include (/home/pkus/sources/hermes-dev/hermes2d/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})


