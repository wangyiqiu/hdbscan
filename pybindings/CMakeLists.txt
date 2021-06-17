project(pyhdbscan)

add_subdirectory(pybind11)

pybind11_add_module(${PROJECT_NAME} pyhdbscan.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE hdbscanLib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
