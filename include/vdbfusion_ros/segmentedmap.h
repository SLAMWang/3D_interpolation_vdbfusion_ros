
#include "utility.h"
#include <pcl/features/normal_3d.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/approximate_voxel_grid.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <Eigen/QR>
using namespace IOU;
struct  Seg_{
  pcl::PointCloud<PointType>::Ptr seg_cloud_;
   Eigen::Vector3f merged_center_;
   Eigen::Matrix3f merged_cov_;
   Eigen::Vector3f last_center_;
   Eigen::Matrix3f last_cov_;
   int obs_times_;
   float dist_;
   
    Vertexes _vert;
};


class SegmentedMap{
 public:

  SegmentedMap(){
    segmented_map_.clear();
    viewpointset.clear();
    frame_index_ = 0;
    global_map_.reset(new pcl::PointCloud<PointType>());
    global_map_weighted_.reset(new pcl::PointCloud<PointType>());
    current_weighted_cloud_.reset(new pcl::PointCloud<PointType>());
    //last_weighted_cloud_.reset(new pcl::PointCloud<PointType>());
    lidar_scans_.clear();
    lidar_scans_seg_.clear();
    scan_odoms_.clear();
    full_lidar_scans_.clear();
    full_cloud_set_.clear();
    kdtree_.reset((new pcl::search::KdTree<PointType>()));
    global_map_pure_.reset(new pcl::PointCloud<PointType>());
    global_ground_map_.reset(new pcl::PointCloud<PointType>());
    lidar_ground_scans_.clear();
    global_ground_cloud_norm_.clear();
      };

  ~SegmentedMap() = default;

  void cloudHandler(pcl::PointCloud<PointType>::Ptr);

 public:
   float computeSegdistance(pcl::PointCloud<PointType>::Ptr cloudtarget,pcl::PointCloud<PointType>::Ptr cloudsource);
   bool SegmentedFrame(pcl::PointCloud<PointType>::Ptr cloud); // Get the segments, 
   void Mapmerge(); //merge the segments into the map
   void Saveviewpoint(PointType p);
   void Mapvisualize();
   bool esti_plane(Eigen::Vector4f &pca_result, std::vector<PointType> points);
   void Mapaggregate(pcl::PointCloud<PointType>::Ptr cloud, pcl::PointCloud<PointType>::Ptr ground_cloud,
   pcl::PointCloud<PointType>::Ptr full_cloud, Eigen::Matrix4f odom);
      void Mapaggregate(pcl::PointCloud<PointType>::Ptr cloud, pcl::PointCloud<PointType>::Ptr ground_cloud, Eigen::Matrix4f odom);
   void Mapaggregate(pcl::PointCloud<PointType>::Ptr cloud,  Eigen::Matrix4f odom);
   void Mapweighted();
   Eigen::MatrixXd  makeScancontext( pcl::PointCloud<PointType> & _scan_down );
   Eigen::MatrixXd  makeMeanScancontext( pcl::PointCloud<PointType> & _scan_down );
   std::vector<std::vector<float>>  detectObs(Eigen::MatrixXd scan_sc, Eigen::MatrixXd map_sc);
   void Getcol(PointType pt, int& row_ind, int& col_ind);
   
  //void publishClouds();
  std::vector<Seg_> segmented_map_;
  std::vector<pcl::PointCloud<PointType>::Ptr> segmented_frame_;
  std::vector<PointType> viewpointset; //each viewpoint xyz of frame
  int frame_index_;
  pcl::PointCloud<PointType>::Ptr global_map_;
  pcl::PointCloud<PointType>::Ptr global_ground_map_;
  pcl::PointCloud<PointType>::Ptr global_map_weighted_;
  pcl::PointCloud<PointType>::Ptr global_map_pure_;
  pcl::PointCloud<PointType>::Ptr current_weighted_cloud_;

  std::vector<Eigen::Vector3d> global_ground_cloud_norm_;

  //std::vector<pcl::PointCloud<PointType>::Ptr> last_scan_quene_;
  pcl::search::KdTree<PointType>::Ptr kdtree_;
  std::vector<pcl::PointCloud<PointType>> lidar_scans_;
  std::vector<pcl::PointCloud<PointType>> lidar_ground_scans_;
  std::vector<pcl::PointCloud<PointType>> full_lidar_scans_;
  std::vector<pcl::PointCloud<PointType>> lidar_scans_seg_;
  std::vector<pcl::PointCloud<PointType>> full_cloud_set_;
  std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> scan_odoms_;

  const double LIDAR_HEIGHT = 1.3; // lidar height : add this for simply directly using lidar scan in the lidar local coord (not robot base coord) / if you use robot-coord-transformed lidar scans, just set this as 0.

    const int    PC_NUM_RING = 20; // 20 in the original paper (IROS 18)
    const int    PC_NUM_SECTOR = 60; // 60 in the original paper (IROS 18)
    const double PC_MAX_RADIUS = 80.0; // 80 meter max in the original paper (IROS 18)
    const double PC_UNIT_SECTORANGLE = 360.0 / double(PC_NUM_SECTOR);
    const double PC_UNIT_RINGGAP = PC_MAX_RADIUS / double(PC_NUM_RING);

    // tree
    const int    NUM_EXCLUDE_RECENT = 50; // simply just keyframe gap, but node position distance-based exclusion is ok. 
    const int    NUM_CANDIDATES_FROM_TREE = 10; // 10 is enough. (refer the IROS 18 paper)

    // loop thres
    const double SEARCH_RATIO = 0.1; // for fast comparison, no Brute-force, but search 10 % is okay. // not was in the original conf paper, but improved ver.
    const double SC_DIST_THRES = 0.13; // empirically 0.1-0.2 is fine (rare false-alarms) for 20x60 polar context (but for 0.15 <, DCS or ICP fit score check (e.g., in LeGO-LOAM) should be required for robustness)
    // const double SC_DIST_THRES = 0.5; // 0.4-0.6 is good choice for using with robust kernel (e.g., Cauchy, DCS) + icp fitness threshold / if not, recommend 0.1-0.15

    // config 
    const int    TREE_MAKING_PERIOD_ = 50; // i.e., remaking tree frequency, to avoid non-mandatory every remaking, to save time cost / if you want to find a very recent revisits use small value of it (it is enough fast ~ 5-50ms wrt N.).
    int          tree_making_period_conter = 0;
    Eigen::MatrixXd depth_sc_;

};

