cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(icp_covariance)

SET(CMAKE_BUILD_TYPE RelWithDebInfo)

find_package(PCL 1.7.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})



add_executable(cov_point_to_point src/point_to_point/loading_inputs_point_to_point.cpp )
target_link_libraries(cov_point_to_point ${PCL_LIBRARIES} )

add_executable(cov_point_to_plane src/point_to_plane/loading_inputs_point_to_plane.cpp )
target_link_libraries(cov_point_to_plane ${PCL_LIBRARIES} )


