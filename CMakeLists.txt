cmake_minimum_required(VERSION 3.12)

project(mc3d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MC3D_SOURCES
    boundary.cpp
    cell.cpp
    cell_claster.cpp
    main.cpp
    p_elem.cpp
    particle.cpp
    particle_array.cpp
    pereodic_boundary.cpp
    point.cpp
    var.cpp
)

add_executable(mc3d ${MC3D_SOURCES})

target_include_directories(mc3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

