cmake_minimum_required(VERSION 2.8)

project(resection)

find_package(Eigen3 REQUIRED)

include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}/src
  )

set(srcs
  ${PROJECT_SOURCE_DIR}/src/photogrammetry.cpp
  )

add_executable(${PROJECT_NAME} ${srcs} "src/main.cpp")

target_link_libraries(${PROJECT_NAME}
  ${EIGEN3_LIBS}
  )
