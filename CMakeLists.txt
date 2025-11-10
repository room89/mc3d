cmake_minimum_required(VERSION 3.16)

project(MC3dSolver VERSION 0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
