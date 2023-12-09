# MIT License
#
# # Copyright (c) 2022 Saurabh Gupta, Ignacio Vizzo, Cyrill Stachniss, University of Bonn
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.0.2)
project(vdbfusion_ros)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# System dependencies are found with CMake's conventions
include(GNUInstallDirs)
list(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/OpenVDB")
find_package(OpenVDB REQUIRED)
find_package(VDBFusion REQUIRED)

link_libraries("/usr/lib/x86_64-linux-gnu/libpcl_features.so")

find_package(
  catkin REQUIRED
  COMPONENTS rospy
             roscpp
             tf2_ros
             pcl_ros
             tf2_sensor_msgs
             std_msgs
             geometry_msgs
             sensor_msgs
             message_generation)

add_service_files(FILES save_vdb_volume.srv)

generate_messages()

find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})


message(STATUS "PCL_INCLUDE_DIRS:" ${PCL_INCLUDE_DIRS})
message(STATUS "PCL_LIBRARY_DIRS:" ${PCL_LIBRARY_DIRS})
message(STATUS "PCL_DEFINITIONS:" ${PCL_DEFINITIONS})


catkin_package(
  INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/include/vdbfusion_ros/
  LIBRARIES
  CATKIN_DEPENDS
  roscpp
  rospy
  tf2_ros
  tf2_sensor_msgs
  std_msgs
  geometry_msgs
  sensor_msgs
  message_runtime)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/vdbfusion_ros ${PCL_INCLUDE_DIRS})


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

include(3rdparty/find_dependencies.cmake)
