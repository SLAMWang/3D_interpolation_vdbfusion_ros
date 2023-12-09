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

add_library(transforms STATIC Transform.cpp imageProjection.cpp segmentedmap.cpp iou.cpp)
target_link_libraries(transforms PUBLIC
  Sophus::Sophus
)
target_include_directories(transforms PRIVATE
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)

SET(CMAKE_BUILD_TYPE Debug)

#add_executable(${PROJECT_NAME}_node VDBVolume_ros_raw.cpp) #for raw vdbfusion
add_executable(${PROJECT_NAME}_node VDBVolume_ros_ustc_final.cpp)
target_link_libraries(${PROJECT_NAME}_node PUBLIC
  ${catkin_LIBRARIES}
  VDBFusion::vdbfusion
  igl::core
  transforms
  ${OpenCV_LIBS}
)
target_include_directories(${PROJECT_NAME}_node PUBLIC
  ${CMAKE_SOURCE_DIR}/include/vdbfusion_ros
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)
