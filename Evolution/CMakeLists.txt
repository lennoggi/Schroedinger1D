cmake_minimum_required(VERSION 3.10)

project(DFT VERSION 1.0.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXE "Schroedinger1D_exe")

add_executable(${EXE}
    Main.cc
    Potentials.cc
    Wave_functions.cc
)

find_package(HDF5 REQUIRED COMPONENTS C)
target_link_libraries(${EXE} PRIVATE HDF5::HDF5)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

install(TARGETS ${EXE}
        RUNTIME DESTINATION bin
)
