// MIT License
//
// # Copyright (c) 2022 Saurabh Gupta, Ignacio Vizzo, Cyrill Stachniss, University of Bonn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>

#include "Transform.hpp"
#include "imageProjection.h"
#include "segmentedmap.h"
#include "vdbfusion/VDBVolume.h"
#include "vdbfusion_ros/save_vdb_volume.h"
#include<boost/bind/bind.hpp>

#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl_ros/point_cloud.h>
#include<pcl/common/transforms.h>
#include<pcl/kdtree/flann.h>
#include<pcl/surface/mls.h>
#include <nav_msgs/Odometry.h>
#include <pcl/filters/random_sample.h>
#include <pcl/filters/impl/random_sample.hpp>
    
#include <pcl/filters/passthrough.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/approximate_voxel_grid.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/statistical_outlier_removal.h>

#include "voxel_map_util.hpp"
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/approximate_voxel_grid.h>
#include<pcl/features/normal_3d.h>

#include <pcl/io/pcd_io.h> 


#include <chrono>
    
//#include "ikd_Tree.hpp"
#define BOOST_BIND_NO_PLACEHOLDERS
#include "ros_time_hash.hpp"

namespace vdbfusion {
class TicToc
{//from A-LOAM, thanks
public:
    TicToc()
    {
        tic();
    }
    explicit TicToc(std::string  _name):name(std::move(_name))
    {
        tic();
    }
    void tic()
    {
        start = std::chrono::system_clock::now();
    }

    double toc()
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count() * 1000;
    }
    void toc_print_ms(){
        std::cout<<name<<":"<< toc() << "ms"<< std::endl;
    }
    void toc_print_us(){
        std::cout<<name<<":"<< toc()*1000 << "us"<< std::endl;
    }

private:
    std::string name{"not named timer"};
    std::chrono::time_point<std::chrono::system_clock> start, end;
};




 typedef pcl::PointXYZI PointT;

class VDBVolumeNode {
    typedef message_filters::sync_policies::ApproximateTime<nav_msgs::Odometry, sensor_msgs::PointCloud2> ApproxSyncPolicy;
    
public:
    VDBVolumeNode();
    std::string kittipose_path; 
    std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>> transpose;
    std::vector<double>  Time_stamp;
    void load_pose_kitti();
    void Integrate();
    Eigen::Matrix4d LiDAT_2_RTK;
    void ComputeChamferDistance(Eigen::MatrixXd V);
    void ComputeTSDFError();





    //VoxelMap
    float max_voxel_size = 0.5;
      int max_layer = 1;
      std::vector<int> layer_size = {5};
      int max_points_size = 100;
      int max_cov_points_size =100;
       float planer_threshold = 0.05;
      std::unordered_map<VOXEL_LOC, PointType> voxel_map_; 
    std::vector<std::pair<PointType,PointType>> augmented_rays_;

private:
    VDBVolume InitVDBVolume();
    void projectPointCloud(pcl::PointCloud<PointT>::Ptr current_cloud,bool update); 
    void initPCLCloud(bool update);
     void checkOcclusion(pcl::PointCloud<PointT>::Ptr cloud);
     size_t computeRangeMapIndex(PointT thisPoint);
    pcl::PointCloud<PointT>::Ptr TriUpSample(PointT A, PointT B, PointT C, float step);
    void UpsamplingPointCloud(pcl::PointCloud<PointT>::Ptr cloud, 
pcl::PointCloud<PointT>::Ptr denser_cloud);
void UpsamplingPointCloudPlane(pcl::PointCloud<PointT>::Ptr cloud, 
pcl::PointCloud<PointT>::Ptr denser_cloud);
    void Integrate(const nav_msgs::OdometryConstPtr& odom_msg, 
    const sensor_msgs::PointCloud2::ConstPtr& cloud_msg, const sensor_msgs::PointCloud2::ConstPtr& submap_msg);
    bool saveVDBVolume(vdbfusion_ros::save_vdb_volume::Request& path,
                       vdbfusion_ros::save_vdb_volume::Response& response);
    pcl::PointCloud<PointT>::Ptr SampleGlobalMap(Eigen::MatrixXd V, Eigen::MatrixXi F);
   // void SampleGlobalMap(Eigen::MatrixXd V, Eigen::MatrixXi F);


    void getSegmentedMap(pcl::PointCloud<PointT>::Ptr cloud);
    void getLidARScans(const sensor_msgs::PointCloud2::ConstPtr& cloud_msg);
    void getLidARPoseScans(const nav_msgs::OdometryConstPtr& odom_msg, const sensor_msgs::PointCloud2::ConstPtr& cloud_msg);
    void BuildVoxelMap();
    float computeMapdistance(pcl::PointCloud<PointT>::Ptr global_map);
    void MapIntegrate();
    void MapIntegrateGlobal();
    void checkDominateCoor(Eigen::Vector3d normvec, Eigen::Vector3d& a,Eigen::Vector3d& b,Eigen::Vector3d& c);
    void AgumentedTSDF(PointType pt,PointType pt0,PointType pt1, std::vector<Eigen::Vector3d>& point_set,
    std::vector<Eigen::Vector3d>& origin_set,Eigen::Vector3d origin );
    void testbaselines();
pcl::search::KdTree<PointT>::Ptr tree_;
pcl::search::KdTree<PointT>::Ptr tree0_;
int interval_ = 3;
bool test_ours_ = true;
bool angle_filter_ground_ = false;
bool angle_filter_ = true;
bool test_vdb = false;
float res_ = 0.2;
private:
    ros::NodeHandle nh_;
    ros::Subscriber sub_;
    ros::ServiceServer srv_;
    Transform tf_;
    ros::Duration timestamp_tolerance_;
    ros::NodeHandle mt_nh;  

    std::unique_ptr<message_filters::Subscriber<nav_msgs::Odometry>> odom_sub;
    std::unique_ptr<message_filters::Subscriber<sensor_msgs::PointCloud2>> cloud_sub;
    std::unique_ptr<message_filters::Subscriber<sensor_msgs::PointCloud2>> submap_sub;
    std::unique_ptr<message_filters::Synchronizer<ApproxSyncPolicy>> sync;
    pcl::PointCloud<PointT>::Ptr laserCloudIn;
    pcl::PointCloud<PointT>::Ptr laserCloudDyn;
    pcl::PointCloud<PointT>::Ptr laserCloudRangeMap;
    Eigen::Isometry3d L_2_I;
    Eigen::Isometry3d last_odom_;
    Eigen::Isometry3d last_keyframe_;
    Eigen::Isometry3d keyframe_;
    pcl::PointCloud<PointT>::Ptr global_map_;
    pcl::PointCloud<PointT>::Ptr Test_global_map_;
    pcl::PointCloud<PointT>::Ptr local_map_;
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> local_rangemap_;
    std::vector<float> spave_curve_dist;
    //KD_TREE<PointT> ikdtree(0.3,0.6,0.2);
    //KD_TREE<PointT> ikd_Tree;
    std::vector<pcl::PointCloud<PointT>::Ptr> keyframe_queue;
    ImageProjection ImPro_;
    SegmentedMap SM_;
    std::vector<pcl::PointCloud<PointType>> lidarscans_;
    std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> scanodoms_;
    
private:
    VDBVolume vdb_volume_;
      int N_SCAN = 16;
     int Horizon_SCAN = 1800;
    // PointCloud Processing
    bool preprocess_;
    bool apply_pose_;
    float min_range_;
    float max_range_;
    int cloud_integrate_num_=0;

    // Triangle Mesh Extraction
    bool fill_holes_;
    float min_weight_;
};
}  // namespace vdbfusion
