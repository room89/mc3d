cmake_minimum_required(VERSION 3.16)

project(MC3dSolver VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

file(GLOB_RECURSE MC3D_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME}
  PRIVATE
    ${MC3D_SOURCES}
)

target_include_directories(${PROJECT_NAME}
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)

if(BUILD_TESTING)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_executable(mc3d_tests
    tests/point_test.cpp
    tests/collision_test.cpp
    tests/cell_cluster_test.cpp
    tests/solver_setup_test.cpp
    tests/conservation_test.cpp
    tests/body_penetration_test.cpp
    src/boundary.cpp
    src/cell_cluster.cpp
    src/cell.cpp
    src/exception.cpp
    src/free_boundary.cpp
    src/geometry.cpp
    src/giper_free_boundary.cpp
    src/inner_boundary.cpp
    src/mirror_boundary.cpp
    src/particle.cpp
    src/pereodic_boundary.cpp
    src/poligon.cpp
    src/point.cpp
    src/solver_setup.cpp
    src/utils/logger.cpp
    src/utils/stopwatch.cpp
    src/utils/threadpool.cpp
    src/utils/utils.cpp
  )

  target_include_directories(mc3d_tests
    PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}/src"
  )

  target_link_libraries(mc3d_tests
    PRIVATE
      GTest::gtest_main
  )

  add_test(NAME mc3d_tests COMMAND mc3d_tests)
endif()
