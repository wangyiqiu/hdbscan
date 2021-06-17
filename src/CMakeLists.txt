project(hdbscan)

add_library(${PROJECT_NAME}Lib
  hdbscan.cpp
  dendrogram.cpp
  )

add_executable(${PROJECT_NAME}
  hdbscanTime.cpp
  hdbscan.cpp
  dendrogram.cpp
  )

# target_link_libraries(${PROJECT_NAME} dendrogram)

#target_link_libraries(${PROJECT_NAME}Lib commonlib dendrogram)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)

target_include_directories(${PROJECT_NAME}Lib PRIVATE ../include)

# target_include_directories(${PROJECT_NAME}Lib PRIVATE ../include)
