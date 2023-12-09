Hi, this is the source code 1.0 of our submitted paper to TIM: Range Map Interpolation based 3D LiDAR Truncated Signed Distance Fields Mapping in Outdoor Environments.

We build our method on vdbfusion, which is an excellent baseline.

Here are some tips to use our code.

1. if you want to know our algorithm, refer to VDBVolume_ros_kitti_maicity_new_final.cpp and VDBVolume_ros_ustc_final.cpp. 
2. If you want to run the code, there exist some absolute paths in our code, you should check and modify them according to yours.
3. The steps to build the environment is the same to vdbfusion.
4. Our code can be implemented on KITTI, Maicity, and our collected VLP-16 dataset.
5. 'roslaunch vdbfusion_ros vdbfusion.launch config_file_name:=KITTI.yaml' is used to collect the lidar scans, then you should run "rosservice call /save_vdb_volume "path: '$Your saved map path$'", and the dense mapping process will be activated.
"

Any problem, please contact me via wangjk@ustc.edu.cn.

Notice that our paper is still under review. 


The following is the readme of vdbfusion. 

# ROS-VDBFusion: Flexible and Efficient TSDF Integration

A ROS C++ wrapper to the [vdbfusion](https://github.com/PRBonn/vdbfusion) library for flexible and efficient TSDF Integration

## Installation

### OpenVDB

```sh
# Install OpenVDB dependencies
sudo apt-get update && apt-get install --no-install-recommends -y \
    libblosc-dev \
    libboost-iostreams-dev \
    libboost-system-dev \
    libboost-system-dev \
    libeigen3-dev

# Install OpenVDB from source
git clone --depth 1 https://github.com/nachovizzo/openvdb.git -b nacho/vdbfusion && cd openvdb
mkdir build && cd build && cmake  -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DUSE_ZLIB=OFF .. &&  sudo make -j$(nproc) all install
```

### VDBfusion

```sh
git clone --depth 1 https://github.com/PRBonn/vdbfusion.git && cd vdbfusion
mkdir build && cd build && cmake .. &&  sudo make -j$(nproc) all install
```

### ROS

For now only [ROS Noetic](http://wiki.ros.org/noetic) is supported.  

As an extra dependency install: `sudo apt install ros-noetic-tf2-sensor-msgs`.
On your catkin workspace, clone the vdbfusion_ros:

```sh
cd catkin_ws/src/
git clone https://github.com/PRBonn/vdbfusion_ros.git
catkin build
```

## Usage

### Configuration

1. Create a config file compatible with your dataset and desired tsdf integration parameters using this [template configuration](/config/template.yaml)
2. The PointClouds should be a `sensor_msgs/PointCloud2` message published on a custom topic name which needs to be specified in the config file
3. The Transforms should be either a `tf2_msgs/TFMessage` or a `geometry_msgs/TransformStamped`. See the [template configuration](config/template.yaml) for more details
4. The data can be published either through a rosbag file or directly from another ros node

### Launch

```sh
roslaunch vdbfusion_ros vdbfusion.launch config_file_name:=<insert config file name here> path_to_rosbag_file:=<insert path to rosbag file here>
```

### Save the VDB Grid and Extract Triangle Mesh

```sh
rosservice call /save_vdb_volume "path: '<insert filename and path to save the volume and mesh>'"    
```

## Dataset Examples

Download the dataset rosbag files from the respective links

### [TU Munich RGB-D SLAM Dataset and Benchmark - FR1DESK2](https://vision.in.tum.de/data/datasets/rgbd-dataset)

- Use the sample [config file](config/FR2Desk2.yaml) provided for this dataset

### [ETH Zurich ASL: Cow and Lady RGBD Dataset](https://projects.asl.ethz.ch/datasets/doku.php?id=iros2017)

- Use the sample [config file](config/CowAndLady.yaml) provided for this dataset

### [KITTI Dataset](http://www.cvlibs.net/datasets/kitti/raw_data.php)

- Convert the dataset into a rosbag file using [kitti2bag](https://github.com/tomas789/kitti2bag)
- Use the sample [config file](config/KITTI.yaml) provided for this dataset

Run the [launch](README.md#launch) command providing config file and rosbag path corresponding to the dataset.

Use the [rosservice](README.md#save-the-vdb-grid-and-extract-triangle-mesh) to save the VDB volume and mesh

## Citation

If you use this library for any academic work, please cite the original [paper](https://www.ipb.uni-bonn.de/wp-content/papercite-data/pdf/vizzo2022sensors.pdf).

```bibtex
@article{vizzo2022sensors,
  author         = {Vizzo, Ignacio and Guadagnino, Tiziano and Behley, Jens and Stachniss, Cyrill},
  title          = {VDBFusion: Flexible and Efficient TSDF Integration of Range Sensor Data},
  journal        = {Sensors},
  volume         = {22},
  year           = {2022},
  number         = {3},
  article-number = {1296},
  url            = {https://www.mdpi.com/1424-8220/22/3/1296},
  issn           = {1424-8220},
  doi            = {10.3390/s22031296}
}
```

## Others ROS wrappers around VDBFusion

- [vdbfusion_mapping](https://github.com/Kin-Zhang/vdbfusion_mapping): originated to solve [transformation issues](https://github.com/PRBonn/vdbfusion_ros/issues/2) and with color support ![image](https://user-images.githubusercontent.com/35365764/200626528-a657a0e6-2fca-48d7-8b34-d8619b6f33e8.png)
- [vdbfusion_mapping_docker](https://github.com/nachovizzo/vdbfusion_mapping). Same as above, but running in a dockerized enviornment

