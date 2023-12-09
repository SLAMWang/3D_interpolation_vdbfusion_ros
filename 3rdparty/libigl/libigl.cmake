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

include(ExternalProject)
ExternalProject_Add(
  external_libigl
  PREFIX libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG v2.3.0
  GIT_SHALLOW ON
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  UPDATE_COMMAND "")

# Simulate importing libigl
ExternalProject_Get_Property(external_libigl SOURCE_DIR)
add_library(libiglHelper INTERFACE)
add_dependencies(libiglHelper external_libigl)
target_include_directories(libiglHelper INTERFACE ${SOURCE_DIR}/include)
add_library(igl::core ALIAS libiglHelper)
