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

FROM osrf/ros:noetic-desktop-full
LABEL maintainer="Ignacio Vizzo <ivizzo@uni-bonn.de>"

# Install utilities
RUN apt-get update && apt-get install --no-install-recommends -y \
    git \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python3 utils
RUN python3 -m pip install --no-cache  catkin-tools

# Install extra ROS dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    ros-${ROS_DISTRO}-tf2-sensor-msgs \
    && rm -rf /var/lib/apt/lists/*

# Install OpenVDB dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    libblosc-dev \
    libboost-iostreams-dev \
    libboost-system-dev \
    libboost-system-dev \
    libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

# Install OpenVDB from source
RUN git clone --depth 1 https://github.com/nachovizzo/openvdb.git -b nacho/vdbfusion \
    && cd openvdb \
    && mkdir build && cd build \
    && cmake  -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DUSE_ZLIB=OFF .. \
    && make -j$(nproc) all install \
    && cd / \
    && rm -rf /openvdb

# Install VDBFusion from source
RUN git clone --depth 1 https://github.com/PRBonn/vdbfusion.git \
    && cd vdbfusion \
    && mkdir build && cd build \
    && cmake .. \
    && make -j$(nproc) all install \
    && cd / \
    && rm -rf /vdbfusion

# Add user to share files and folder without root permissions
ENV GROUP_ID=1000
ENV USER_ID=1000
RUN addgroup --gid $GROUP_ID user && adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user
USER user

# Build ROS workspace, bash needs the "-l" to load the /etc/profile
WORKDIR /home/user/ros_ws

ENTRYPOINT ["/ros_entrypoint.sh"]
CMD ["bash", "--login"]
