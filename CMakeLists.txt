cmake_minimum_required(VERSION 3.16)

project(MC3dSolver VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

add_library(mc3d_core STATIC
  src/boundary.cpp
  src/collision_model.cpp
  src/particle_sampler.cpp
  src/cell.cpp
  src/cell_cluster.cpp
  src/cellnet.cpp
  src/config.cpp
  src/run_manifest.cpp
  src/exception.cpp
  src/free_boundary.cpp
  src/geometry.cpp
  src/giper_free_boundary.cpp
  src/inner_boundary.cpp
  src/inner_cell.cpp
  src/mirror_boundary.cpp
  src/particle.cpp
  src/pereodic_boundary.cpp
  src/point.cpp
  src/poligon.cpp
  src/simulation_runner.cpp
  src/solver_setup.cpp
  src/utils/logger.cpp
  src/utils/stopwatch.cpp
  src/utils/threadpool.cpp
  src/utils/utils.cpp
  src/var.cpp
  src/warning.cpp
)
target_include_directories(mc3d_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(mc3d_core PUBLIC Threads::Threads)
add_custom_target(mc3d_build_revision
  COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_SOURCE_DIR} -DBINARY_DIR=${CMAKE_BINARY_DIR}
    -P ${CMAKE_SOURCE_DIR}/tools/write_build_revision.cmake
  BYPRODUCTS ${CMAKE_BINARY_DIR}/generated/build_revision.h
  VERBATIM)
add_dependencies(mc3d_core mc3d_build_revision)
configure_file(src/build_info.h.in generated/build_info.h @ONLY)
target_include_directories(mc3d_core PRIVATE "${CMAKE_BINARY_DIR}/generated")
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE mc3d_core)

if(BUILD_TESTING)
  include(FetchContent)
  include(GoogleTest)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_executable(mc3d_tests
    tests/geometry_queries_test.cpp
    tests/random_test.cpp
    tests/point_test.cpp
    tests/value_types_test.cpp
    tests/collision_scattering_test.cpp
    tests/collision_frequency_test.cpp
    tests/vhs_test.cpp
    tests/collision_test.cpp
    tests/collision_normalization_test.cpp
    tests/cell_cluster_test.cpp
    tests/solver_setup_test.cpp
    tests/conservation_test.cpp
    tests/surface_heating_test.cpp
    tests/body_penetration_test.cpp
    tests/open_boundary_test.cpp
    tests/hyperfree_test.cpp
  )

  target_include_directories(mc3d_tests
    PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}/src"
  )

  target_link_libraries(mc3d_tests
    PRIVATE
      mc3d_core
      GTest::gtest_main
  )

  gtest_discover_tests(mc3d_tests)
endif()

if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_test(NAME RunManifestIntegration COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/tests/run_manifest_integration.py $<TARGET_FILE:${PROJECT_NAME}>)
  set_tests_properties(RunManifestIntegration PROPERTIES TIMEOUT 120)
endif()
if(BUILD_TESTING)
  add_test(NAME ConfigExamplesIntegration COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/tests/config_examples_integration.py $<TARGET_FILE:${PROJECT_NAME}>)
  set_tests_properties(ConfigExamplesIntegration PROPERTIES TIMEOUT 120)
endif()
