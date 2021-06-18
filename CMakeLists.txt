cmake_minimum_required(VERSION 3.10)

project(parallel-hdbscan)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-DHOMEGROWN -pthread -O3 -DNDEBUG -fPIC ${EXTRA_FLAGS}")

add_subdirectory(src)

add_subdirectory(pybindings)
