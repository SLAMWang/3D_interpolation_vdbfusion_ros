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

#include "VDBVolume_ros.hpp"
#include <ctime>
#include <geometry_msgs/TransformStamped.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/point_cloud_conversion.h>
#include <tf/transform_listener.h>
#include <tf2_sensor_msgs/tf2_sensor_msgs.h>
#include <boost/filesystem.hpp>

#include <Eigen/Core>
#include <vector>
#include <pcl/io/ply_io.h>
#include "igl/write_triangle_mesh.h"
#include "openvdb/openvdb.h"


namespace {
          typedef pcl::PointXYZI PointT;

std::vector<Eigen::Vector3d> pcl2SensorMsgToEigen(pcl::PointCloud<PointT>::Ptr pcl) {
    std::vector<Eigen::Vector3d> points;
    points.reserve(pcl->size());
    std::for_each(pcl->points.begin(), pcl->points.end(), [&](const auto& point) {
        points.emplace_back(Eigen::Vector3d(point.x, point.y, point.z));
    });
    return points;
}

void PreProcessCloud(std::vector<Eigen::Vector3d>& points, float min_range, float max_range) {
    points.erase(
        std::remove_if(points.begin(), points.end(), [&](auto p) { return p.norm() > max_range; }),
        points.end());
    points.erase(
        std::remove_if(points.begin(), points.end(), [&](auto p) { return p.norm() < min_range; }),
        points.end());
}
}  // namespace

vdbfusion::VDBVolume vdbfusion::VDBVolumeNode::InitVDBVolume() {
    float voxel_size;
    float sdf_trunc;
    bool space_carving;

    nh_.getParam("/voxel_size", voxel_size);
    nh_.getParam("/sdf_trunc", sdf_trunc);
    nh_.getParam("/space_carving", space_carving);

    VDBVolume vdb_volume(voxel_size, sdf_trunc, space_carving);
    local_map_.reset(new pcl::PointCloud<PointT>());
    return vdb_volume;
}

vdbfusion::VDBVolumeNode::VDBVolumeNode() : vdb_volume_(InitVDBVolume()), tf_(nh_) {
    openvdb::initialize();
    std::string pcl_topic;
    std::cout<<"Initialize"<<std::endl;
    nh_.getParam("/pcl_topic", pcl_topic);
    nh_.getParam("/preprocess", preprocess_);
    nh_.getParam("/apply_pose", apply_pose_);
    nh_.getParam("/min_range", min_range_);
    nh_.getParam("/max_range", max_range_);

    nh_.getParam("/fill_holes", fill_holes_);
    nh_.getParam("/min_weight", min_weight_);

    int32_t tol;
    nh_.getParam("/timestamp_tolerance_ns", tol);
    timestamp_tolerance_ = ros::Duration(0, tol);
    std::cout<<"Initialize"<<std::endl;
    const int queue_size = 500;

    L_2_I= Eigen::Isometry3d::Identity();
    Eigen::Matrix3d Lidar_R_wrt_IMU;
    Eigen::Vector3d Lidar_T_wrt_IMU;
    Lidar_R_wrt_IMU<<0.999371, -0.0353897, 0.00226962,
                -0.0354216, -0.999243, 0.0160615,
                 0.0016995, -0.0161318, -0.999868;
    Lidar_T_wrt_IMU<<-0.0156031,  -0.028666, -0.661778;
    L_2_I.translate(Lidar_T_wrt_IMU);
    L_2_I.rotate(Lidar_R_wrt_IMU);
    
    std::cout<<"Initialize"<<std::endl;
    odom_sub.reset(new message_filters::Subscriber<nav_msgs::Odometry>(nh_, "/Odometry", 256));
    cloud_sub.reset(new message_filters::Subscriber<sensor_msgs::PointCloud2>(nh_, "/velodyne_points", 32));
   
    sub_ = nh_.subscribe("/velodyne_points", 32, &vdbfusion::VDBVolumeNode::getLidARScans, this);
    //sub_ = nh_.subscribe(pcl_topic, queue_size, &vdbfusion::VDBVolumeNode::Integrate, this);
    srv_ = nh_.advertiseService("/save_vdb_volume", &vdbfusion::VDBVolumeNode::saveVDBVolume, this);
   std::cout<<"Initialize"<<std::endl;
    laserCloudRangeMap.reset(new pcl::PointCloud<PointT>());
    laserCloudRangeMap->resize(N_SCAN*Horizon_SCAN);
     PointT  nanPoint;
   nanPoint.x = std::numeric_limits<float>::quiet_NaN();
   nanPoint.y = std::numeric_limits<float>::quiet_NaN();
   nanPoint.z = std::numeric_limits<float>::quiet_NaN();
   nanPoint.intensity = std::numeric_limits<float>::quiet_NaN();
    std::fill(laserCloudRangeMap->points.begin(), laserCloudRangeMap->points.end(), nanPoint);
    ROS_INFO("Use '/save_vdb_volume' service to save the integrated volume");
    laserCloudDyn.reset(new pcl::PointCloud<PointT>());
    transpose.clear();
    Time_stamp.clear();
    LiDAT_2_RTK << 4.276802385584e-04, -9.999672484946e-01, -8.084491683471e-03, -1.198459927713e-02,
            -7.210626507497e-03, 8.081198471645e-03, -9.999413164504e-01, -5.403984729748e-02,
            9.999738645903e-01, 4.859485810390e-04, -7.206933692422e-03, -2.921968648686e-01,
            0, 0, 0, 1;
    tree_.reset(new pcl::search::KdTree<PointT>());
    cloud_integrate_num_ = 0;
    Test_global_map_.reset(new pcl::PointCloud<PointT>());
}

static Eigen::Isometry3d odom2isometry(const nav_msgs::OdometryConstPtr& odom_msg) {
  const auto& orientation = odom_msg->pose.pose.orientation;
  const auto& position = odom_msg->pose.pose.position;

  Eigen::Quaterniond quat;
  quat.w() = orientation.w;
  quat.x() = orientation.x;
  quat.y() = orientation.y;
  quat.z() = orientation.z;

  Eigen::Isometry3d isometry = Eigen::Isometry3d::Identity();
  isometry.linear() = quat.toRotationMatrix();
  isometry.translation() = Eigen::Vector3d(position.x, position.y, position.z);
  return isometry;
}




void vdbfusion::VDBVolumeNode::initPCLCloud(bool localmap_update_){
   
   laserCloudIn.reset(new pcl::PointCloud<PointT>());
    laserCloudIn->resize(N_SCAN*Horizon_SCAN);
    PointT  nanPoint;
   nanPoint.x = std::numeric_limits<float>::quiet_NaN();
   nanPoint.y = std::numeric_limits<float>::quiet_NaN();
   nanPoint.z = std::numeric_limits<float>::quiet_NaN();
   nanPoint.intensity = std::numeric_limits<float>::quiet_NaN();
    std::fill(laserCloudIn->points.begin(), laserCloudIn->points.end(), nanPoint);
    
    
    if(localmap_update_)
    {
        Eigen::Vector3d pt(10000,0,0); //mean cov, num
        local_rangemap_.resize(N_SCAN*Horizon_SCAN);
        std::fill(local_rangemap_.begin(),local_rangemap_.end(),pt);
    }
         
    
}

size_t vdbfusion::VDBVolumeNode::computeRangeMapIndex(PointT thisPoint)
{
                float ang_res_x = 0.2*M_PI/180;
                float ang_res_y = 2.0*M_PI/180;
                float ang_bottom = (15.0+0.1)*M_PI/180;
                float range = sqrt(thisPoint.x * thisPoint.x +
                                    thisPoint.y * thisPoint.y +
                                    thisPoint.z * thisPoint.z);
                 float verticalAngle = std::asin(thisPoint.z / range);
                    //std::atan2(thisPoint.z, sqrt(thisPoint.x * thisPoint.x + thisPoint.y * thisPoint.y));

                int rowIdn = (verticalAngle + ang_bottom) / ang_res_y;
                if (rowIdn < 0 || rowIdn >= N_SCAN) {
                return -1;
                }

                float horizonAngle = std::atan2(thisPoint.x, thisPoint.y);

                int columnIdn = -round((horizonAngle - M_PI_2) / ang_res_x) + Horizon_SCAN * 0.5;

                if (columnIdn >= Horizon_SCAN){
                columnIdn -= Horizon_SCAN;
                }

                if (columnIdn < 0 || columnIdn >= Horizon_SCAN){
                return -1;
                }

                if (range < 0.1){
                return -1;
                }
               size_t index = columnIdn + rowIdn * Horizon_SCAN;
               return index;
}

   void vdbfusion::VDBVolumeNode::projectPointCloud(pcl::PointCloud<PointT>::Ptr current_cloud, bool localmap_update_) {
             initPCLCloud(localmap_update_);
            // range image projection
            // ROS msg -> PCL cloud
            
            
            pcl::PointCloud<PointT>::Ptr laserCloudRawIn(new pcl::PointCloud<PointT>());
            pcl::copyPointCloud(*current_cloud,*laserCloudRawIn);
            for (size_t i = 0; i < laserCloudRawIn->size(); ++i) {
                PointT thisPoint = laserCloudRawIn->points[i];

                float range = sqrt(thisPoint.x * thisPoint.x +
                                thisPoint.y * thisPoint.y +
                                thisPoint.z * thisPoint.z);
                 size_t index = computeRangeMapIndex(thisPoint);
                // find the row and column index in the image for this point
                if(index==-1)
                    continue;
                // the corresponding range of a point is saved as "intensity"
                thisPoint.intensity = range;
                if(isnan(laserCloudIn->points[index].intensity))
                    laserCloudIn->points[index] = thisPoint;
                else
                {
                     if(range< laserCloudIn->points[index].intensity)
                               laserCloudIn->points[index] = thisPoint;
                }
            }
}

float kernal(Eigen::Vector3f X1,Eigen::Vector3f X2)
{
      return exp(-1*pow((X1-X2).norm(),2)/(2*pow(0.5,2)));
}

pcl::PointCloud<PointT>::Ptr vdbfusion::VDBVolumeNode::TriUpSample(PointT A, PointT B, PointT C, float step)
{
     pcl::PointCloud<PointT>::Ptr laserCloud(new pcl::PointCloud<PointT>());
      PointT New_point;

      if(step<0)
      {
           laserCloud->points.push_back(A);
           laserCloud->points.push_back(B);
           laserCloud->points.push_back(C);
           laserCloud->width = 1;
           laserCloud->height = 3;
           return laserCloud;
      }

      bool gaussian_proc = false;
      if(gaussian_proc)
      {
        Eigen::Vector3f a = A.getVector3fMap() - B.getVector3fMap();
     Eigen::Vector3f b = C.getVector3fMap() - B.getVector3fMap();
    Eigen::Vector3f n = a.cross(b);
    n.normalize();
           std::ptrdiff_t ind;
           
           float maxOfV = n.cwiseAbs().maxCoeff(&ind);
           //std::cout<<"n is "<<n.transpose()<<" "<<maxOfV<<" "<<ind<<std::endl;

           Eigen::Vector3f Av = A.getVector3fMap();
           Eigen::Vector3f Bv = B.getVector3fMap();
           Eigen::Vector3f Cv = C.getVector3fMap();

           Eigen::Vector3f Av0 = A.getVector3fMap();
           Eigen::Vector3f Bv0 = B.getVector3fMap();
           Eigen::Vector3f Cv0 = C.getVector3fMap();

           Av0[ind] = 0;Bv0[ind] = 0;Cv0[ind] = 0;
           Eigen::Vector3f Dv0;
           int num = 0;
           //std::cout<<"Av: "<<Av.transpose()<<std::endl;
           //std::cout<<"Bv: "<<Bv.transpose()<<std::endl;
           //std::cout<<"Cv: "<<Cv.transpose()<<std::endl;
            for(float i = 0.1; i < 1; i = i+step)
                {
                    for(float j = 0.1; j < 1; j= j+step)
                    {
                        float k = 1-i-j;
                        if(k<0)
                            continue;
                        Dv0 = i * Av0 + j * Bv0 + k * Cv0;
                        Eigen::Matrix3f K;
                        K<<kernal(Av0,Av0),kernal(Av0,Bv0),kernal(Av0,Cv0),
                                kernal(Bv0,Av0),kernal(Bv0,Bv0),kernal(Bv0,Cv0),
                                kernal(Cv0,Av0),kernal(Cv0,Bv0),kernal(Cv0,Cv0);
                        Eigen::Vector3f K0(kernal(Dv0,Av0),kernal(Dv0,Bv0),kernal(Dv0,Cv0));

                        float y_pred = K0.transpose() * (K+0.1*Eigen::Matrix3f::Identity()).inverse() * Eigen::Vector3f(Av[ind],Bv[ind],Cv[ind]);

                        Dv0[ind] = y_pred;
                        New_point.getVector3fMap() = Dv0;
                        laserCloud->points.push_back(New_point);
                        num++;
                        //std::cout<<"Dv0 "<<Dv0.transpose()<<std::endl;
                    }
                }
           

     laserCloud->width = 1;
     laserCloud->height = num;
     return laserCloud;



      }
     
     bool idv = false;
     if(idv)
     {
         Eigen::Vector3f a = A.getVector3fMap() - B.getVector3fMap();
        Eigen::Vector3f b = C.getVector3fMap() - B.getVector3fMap();
        Eigen::Vector3f n = a.cross(b);
        n.normalize();
           std::ptrdiff_t ind;
           
           float maxOfV = n.cwiseAbs().maxCoeff(&ind);
           //std::cout<<"n is "<<n.transpose()<<" "<<maxOfV<<" "<<ind<<std::endl;

           Eigen::Vector3f Av = A.getVector3fMap();
           Eigen::Vector3f Bv = B.getVector3fMap();
           Eigen::Vector3f Cv = C.getVector3fMap();

           Eigen::Vector3f Av0 = A.getVector3fMap();
           Eigen::Vector3f Bv0 = B.getVector3fMap();
           Eigen::Vector3f Cv0 = C.getVector3fMap();

           Av0[ind] = 0;Bv0[ind] = 0;Cv0[ind] = 0;
           Eigen::Vector3f Dv0;
           int num = 0;
           //std::cout<<"Av: "<<Av.transpose()<<std::endl;
           //std::cout<<"Bv: "<<Bv.transpose()<<std::endl;
           //std::cout<<"Cv: "<<Cv.transpose()<<std::endl;
            for(float i = 0.1; i < 1; i = i+step)
                {
                    for(float j = 0.1; j < 1; j= j+step)
                    {
                        float k = 1-i-j;
                        if(k<0)
                            continue;
                        Dv0 = i * Av0 + j * Bv0 + k * Cv0;
                        Eigen::Matrix3f K;
                        
                        float d0 = (Av0 - Dv0).norm();
                        float d1 = (Bv0 - Dv0).norm();
                        float d2 = (Cv0 - Dv0).norm();

                        float y_pred = (Av[ind]/d0 + Bv[ind]/d1 + Cv[ind]/d2)/(1/d0 + 1/d1 + 1/d2);

                        Dv0[ind] = y_pred;
                        New_point.getVector3fMap() = Dv0;
                        laserCloud->points.push_back(New_point);
                        num++;
                        //std::cout<<"Dv0 "<<Dv0.transpose()<<std::endl;
                    }
                }
           

     laserCloud->width = 1;
     laserCloud->height = num;
     return laserCloud;
     }

     
     
     int num = 0;
     for(float i = 0.1; i < 1; i = i+step)
     {
         for(float j = 0.1; j < 1; j= j+step)
         {
             float k = 1-i-j;
             if(k<0)
                continue;
             
             New_point.x = i * A.x + j * B.x + k * C.x;
             New_point.y = i * A.y + j * B.y + k * C.y;
             New_point.z = i * A.z + j * B.z + k * C.z;
             New_point.intensity = A.intensity;
             laserCloud->points.push_back(New_point);
             num++;
         }
     }
     laserCloud->width = 1;
     laserCloud->height = num;
     return laserCloud;
}

Eigen::Vector3f ComputeMeshNormal(PointT A, PointT B, PointT C)
{
     Eigen::Vector3f a(A.x-B.x,A.y-B.y,A.z-B.z);
     Eigen::Vector3f b(C.x-B.x,C.y-B.y,C.z-B.z);
    Eigen::Vector3f c = a.cross(b);
     return c;
}

void vdbfusion::VDBVolumeNode::UpsamplingPointCloudPlane(pcl::PointCloud<PointT>::Ptr cloud, 
pcl::PointCloud<PointT>::Ptr denser_cloud)
{
    //TODO: firstly, generate rangemap; secondly, obtain each rectangle
    //TODO: thirdly, for each retangle, judge its effectiveness, and then perform interplotation
    *denser_cloud += *cloud;
    projectPointCloud(cloud,false); //obtain the laserCloudIn
    denser_cloud->clear();
    int k = 2;
    for (int i = 1; i < N_SCAN-1; i++)
    {
         for(int j = 1; j < Horizon_SCAN-5; j=j+k)
         {
             size_t ind0 = i * Horizon_SCAN + j;
             size_t ind1 = i * Horizon_SCAN + j+k;
             size_t ind2 = (i+1) * Horizon_SCAN + j;
             size_t ind3 = (i+1) * Horizon_SCAN + j+k;
             
            bool up_added = false;
            bool down_added = false;
            pcl::PointCloud<PointT>::Ptr  samples_cloud0(new pcl::PointCloud<PointT>());
            pcl::PointCloud<PointT>::Ptr  samples_cloud1(new pcl::PointCloud<PointT>());
            samples_cloud0->clear();
            samples_cloud1->clear();
             //std::cout<<"range: "<<laserCloudIn->points[ind0].intensity<<" "<<laserCloudIn->points[ind1].intensity<<" "<<
             //laserCloudIn->points[ind2].intensity<<std::endl;
            if(isnan(laserCloudIn->points[ind0].intensity))
                 continue;
             
             if(isnan(laserCloudIn->points[ind1].intensity))
                  continue;
            
             if(isnan(laserCloudIn->points[ind2].intensity))
                  continue;
             Eigen::Vector3f n1 = ComputeMeshNormal(laserCloudIn->points[ind0],
              laserCloudIn->points[ind1], laserCloudIn->points[ind2]);
              Eigen::Vector3f n2 = ComputeMeshNormal(laserCloudIn->points[ind1],
              laserCloudIn->points[ind2], laserCloudIn->points[ind3]);


              
              Eigen::Vector3f view0 = laserCloudIn->points[ind0].getVector3fMap();
               Eigen::Vector3f view3 = laserCloudIn->points[ind3].getVector3fMap();
               /*
              if(abs(n1.dot( view0 ))/(n1.norm()* view0 .norm())>0.8)
              {
                        samples_cloud0 = TriUpSample(laserCloudIn->points[ind0],
                            laserCloudIn->points[ind1], laserCloudIn->points[ind2],0.2);
                        up_added = true;
              }
              if(abs(n2.dot( view3 ))/(n2.norm()* view3 .norm())>0.8)
              {
                         samples_cloud1 = TriUpSample(laserCloudIn->points[ind1],
                            laserCloudIn->points[ind2], laserCloudIn->points[ind3],0.2);
                        down_added = true;
              }
              if(abs(n1.dot( view0 ))/(n1.norm()* view0 .norm())<0.2)
                   continue;
  */

              if(abs(n1.dot(n2))/(n1.norm()*n2.norm())<0.95)
                   continue;
              if(abs(n1.dot( view0 ))/(n1.norm()* view0 .norm())<0.15)
                   continue;

             if(true )
             {
                 samples_cloud0   = TriUpSample(laserCloudIn->points[ind0],
                 laserCloudIn->points[ind1], laserCloudIn->points[ind2],0.4);
             }
                
//
            
             if(isnan(laserCloudIn->points[ind3].intensity))
                  continue;
            if(true)
                   samples_cloud1 = TriUpSample(laserCloudIn->points[ind1],
              laserCloudIn->points[ind2], laserCloudIn->points[ind3],0.4);
             *denser_cloud +=  *samples_cloud0;
             *denser_cloud +=  *samples_cloud1;
         }
    }
}

void vdbfusion::VDBVolumeNode::UpsamplingPointCloud(pcl::PointCloud<PointT>::Ptr cloud, 
pcl::PointCloud<PointT>::Ptr denser_cloud)
{
    //TODO: firstly, generate rangemap; secondly, obtain each rectangle
    //TODO: thirdly, for each retangle, judge its effectiveness, and then perform interplotation
    *denser_cloud += *cloud;
    projectPointCloud(cloud,false); //obtain the laserCloudIn
    int k = 2;
    for (int i = 0; i < N_SCAN-1; i++)
    {
         for(int j = 1; j < Horizon_SCAN-3; j=j+k)
         {
             size_t ind0 = i * Horizon_SCAN + j;
             size_t ind2 = (i+1) * Horizon_SCAN + j;             
            bool up_added = false;
            bool down_added = false;
            pcl::PointCloud<PointT>::Ptr  samples_cloud0(new pcl::PointCloud<PointT>());
            samples_cloud0->clear();
            if(isnan(laserCloudIn->points[ind0].intensity))
                 continue;
             
            
             if(isnan(laserCloudIn->points[ind2].intensity))
                  continue;
             
             Eigen::Vector3f a = laserCloudIn->points[ind0].getVector3fMap();
             Eigen::Vector3f b = laserCloudIn->points[ind2].getVector3fMap();
             float dist_c = (a-b).norm();
             float dist_a = laserCloudIn->points[ind0].intensity;
             float dist_b = laserCloudIn->points[ind2].intensity;
             float costheta = (pow(dist_c,2)+pow(dist_a,2)-pow(dist_b,2))/(2*dist_c*dist_a);
             if(abs(costheta)>0.35)
                 continue;
             //if(abs(n[2])>0.1 and abs(n[2])<0.9)
             //    continue;
             
             if(laserCloudIn->points[ind0].intensity>50)
                  continue;

              PointT C;
              C.x = 0;
              C.y = 0;
              C.z = 0;
              C.intensity = 0;
                
                 samples_cloud0   = TriUpSample(laserCloudIn->points[ind0],
                 laserCloudIn->points[ind2], C,0.2);

             *denser_cloud +=  *samples_cloud0;
         }
    }
}

void vdbfusion::VDBVolumeNode::checkOcclusion(pcl::PointCloud<PointT>::Ptr cloud)
{
    //first, transform the last rangemap into the current pose
    //recompute the rangemap
    //compare and update the rangemap

    //directly add the global point cloud but the problem is the map cloud will be size increased. 
    laserCloudDyn->clear();
        spave_curve_dist.resize(cloud->size(),1.0);
        int c = 0;
        std::cout<<"cloud size: "<<cloud->size()<<std::endl;
        for(int i = 0; i < cloud->size(); i++)
        {
            PointT thisPoint = cloud->points[i];

            size_t index = computeRangeMapIndex(thisPoint);
            //std::cout<<"index "<<index<<" "<<thisPoint.x<<" "<<thisPoint.y<<" "<<
            //thisPoint.z<<std::endl;
            float range = sqrt(thisPoint.x * thisPoint.x +
                                    thisPoint.y * thisPoint.y +
                                    thisPoint.z * thisPoint.z);
            
            //Eigen::Vector3d mv = local_rangemap_.at(index);
            if(abs(range-laserCloudIn->points[index].intensity)>0.25)
            {
                laserCloudDyn->points.push_back(thisPoint);
                    spave_curve_dist.at(i) = 1;
                    c++;
            }else
                   spave_curve_dist.at(i) = -1;
                 
        }
        std::cout<<"c is "<<c<<std::endl;

}
//TODO 1. select keyframe; 2. aggregate the local map; 3. project the rangemap; 4. select the points that need to be space curving




void vdbfusion::VDBVolumeNode::getLidARScans(const sensor_msgs::PointCloud2::ConstPtr& cloud_msg)
{
     pcl::PointCloud<PointT>::Ptr cloud_raw(new pcl::PointCloud<PointT>());
     pcl::fromROSMsg(*cloud_msg, *cloud_raw);
     lidarscans_.push_back(*cloud_raw);
}



PointT point2pcl(Eigen::VectorXd point)
{
     PointT pt;
     pt.x = point[0];
     pt.y = point[1];
     pt.z = point[2];
     pt.intensity =0;
     return pt;
}



float computeMapdistance(pcl::PointCloud<PointT>::Ptr global_map)
{
    pcl::PointCloud<PointT>::Ptr ref_Map(new pcl::PointCloud<PointT>());
    //std::string ref_map_path = "/home/wjk/Fast-LIO2/src/FAST_LIO/PCD/scans.pcd";
    std::string ref_map_path =  "/home/wjk/catkin_ws_vbdfusion/results/globalmap.pcd";
        if(-1 == pcl::io::loadPCDFile(ref_map_path,*ref_Map)){
        std::cout<<"load pcd file failed. please check it."<<std::endl;
        return -2;
    }

  pcl::search::KdTree<PointT>::Ptr tree_(new pcl::search::KdTree<PointT>());
  tree_->setInputCloud(global_map);

  double fitness_score = 0.0;

  std::vector<int> nn_indices(1);
  std::vector<float> nn_dists(1);

   pcl::PointCloud<PointT>::Ptr cloud_sub(new pcl::PointCloud<PointT>);    //随机下采样点云
    pcl::RandomSample<PointT> rs;    //创建滤波器对象
    rs.setInputCloud(ref_Map);                //设置待滤波点云
    rs.setSample(200000);                    //设置下采样点云的点数
    //rs.setSeed(1);                        //设置随机函数种子点
    rs.filter(*cloud_sub);   

  // For each point in the source dataset
  int nr = 0;
  int point_num = 1;
  for(size_t i = 0; i < cloud_sub->size(); ++i) {
    // Find its nearest neighbor in the target
    tree_->nearestKSearch(cloud_sub->points[i], point_num, nn_indices, nn_dists);
    if( nn_dists[0] <= 0.5) {
      // Add to the fitness score
      fitness_score += pow(nn_dists[0],1);
      nr++;
    }
  }
  if(nr > 0)
    return (fitness_score / nr);
  else
    return (std::numeric_limits<double>::max());
}



void vdbfusion::VDBVolumeNode::ComputeChamferDistance(Eigen::MatrixXd V)
{
    pcl::PointCloud<PointT>::Ptr V_Map(new pcl::PointCloud<PointT>());
 

    
    tree_->setInputCloud(Test_global_map_);
    std::vector<int> nn_indices(1);
    std::vector<float> nn_dists(1);

    std::cout<<"V.rows"<<std::endl;
    std::cout<<V.rows()<<std::endl;
    std::cout<< V.row(0)(0)<<" "<< V.row(0)(1)<<" "<< V.row(0)(2)<<std::endl;
    float chamfer_dist0 = 0;
    V_Map->points.clear();
    for(int i = 0; i < V.rows(); i++)
    {
        PointType pt;
        pt.x = V.row(i)(0);
        pt.y = V.row(i)(1);
        pt.z = V.row(i)(2);
        pt.intensity = 0;
        V_Map->points.push_back(pt);
        tree_->nearestKSearch(pt, 1, nn_indices, nn_dists);
        chamfer_dist0 += nn_dists[0];
    }
    chamfer_dist0 = chamfer_dist0/(2*V.rows());
     
    V_Map->height =  V.rows();
    V_Map->width = 1;
    std::cout<<"chamfer_dist0: "<<chamfer_dist0<<std::endl;
     tree_->setInputCloud(V_Map);
     float chamfer_dist1 = 0;
    for(int i = 0; i < Test_global_map_->size(); i++)
    {
         PointType pt = Test_global_map_->points[i];
         if(isnan(pt.x) or isnan(pt.y) or isnan(pt.z))
            continue;
         tree_->nearestKSearch(pt, 1, nn_indices, nn_dists);
         chamfer_dist1 += nn_dists[0];

    }
    chamfer_dist1 = chamfer_dist1/(2*Test_global_map_->size());
    std::cout<<"chamfer_dist1: "<<chamfer_dist1<<" "<<Test_global_map_->size()<<std::endl;
    std::cout<<"chamfer_dist: "<<chamfer_dist0+chamfer_dist1<<std::endl;
    
}

void vdbfusion::VDBVolumeNode::ComputeTSDFError()
{
        //std::string ref_map_path = "/home/wjk/Fast-LIO2/src/FAST_LIO/PCD/scans.pcd";
        //std::string ref_map_path =  "/home/wjk/catkin_ws_vbdfusion/results/globalmap.pcd";
        /*
        std::string ref_map_path =  "/home/wjk/Dataset/MaiCity/mai_city/gt_models/mai_city_block.ply";
        if(-1 == pcl::io::loadPLYFile(ref_map_path,*ref_Map)){
        std::cout<<"load pcd file failed. please check it."<<std::endl;
            return;
    }*/
    
        
    float tsdf_error = 0;
    for(size_t i = 0; i < Test_global_map_->size(); ++i)
    {
         PointType  pt = Test_global_map_->points[i];
         Eigen::Vector3d point(pt.x,pt.y,pt.z);
         float s = abs(vdb_volume_.GetTSDFValue(point));
         if(s>0.98)
            tsdf_error += 10;
        else
             tsdf_error += s;

    }
    std::cout<<"tsdf_error: "<<(float)tsdf_error/Test_global_map_->size()<<std::endl;
}

void vdbfusion::VDBVolumeNode::load_pose_kitti() {
        ifstream pose(kittipose_path);
        cout << "path: " << kittipose_path << endl;
        std::string line;
        while (getline(pose, line)) {
            if (line.size() > 0) {
                std::stringstream ss(line);
                Eigen::Matrix4d Pose_(Eigen::Matrix4d::Identity());
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 4; ++j) {
                        ss >> Pose_(i, j);
                    }
                transpose.push_back(Pose_);
            }
        }
        pose.close();
        /*
        ofstream tum_pose(sequence+".txt");
        for (int i=0;i<transpose.size();i++){
            Eigen::Quaterniond quaterniond;
            quaterniond=transpose[i].block<3,3>(0,0);
            tum_pose<<Time_stamp[i]<<" "<<transpose[i](0,3)<<" "<<transpose[i](1,3)<<" "<<transpose[i](2,3)<<" "
                    <<quaterniond.x()<<" "<<quaterniond.y()<<" "<<quaterniond.z()<<" "<<quaterniond.w()<<endl;
        }
        tum_pose.close();*/

    }

void vdbfusion::VDBVolumeNode::Integrate()
{
    std::cout<<"lidarscans_ "<<lidarscans_.size()<<" "<<transpose.size()<<std::endl;
    pcl::PointCloud<PointT>::Ptr global_cloud(new pcl::PointCloud<PointT>());
    pcl::PointCloud<PointT>::Ptr global_ground_cloud(new pcl::PointCloud<PointT>());
    pcl::PointCloud<PointT>::Ptr global_full_cloud(new pcl::PointCloud<PointT>());
    Eigen::Matrix4d T0 = transpose.at(0)* LiDAT_2_RTK;
      for(int i = 0; i  <std::min(1000,(int)lidarscans_.size())-2; i=i+2)
      {
        std::cout<<"integrate "<<i<<std::endl;
        ImPro_.cloudHandler(lidarscans_.at(i).makeShared());
        
           global_cloud->clear();
           Eigen::Matrix4d T1 =  transpose.at(i+1)* LiDAT_2_RTK;//// kitti
            //Eigen::Matrix4d T1 =  transpose.at(i+1);
            //std::cout<<"T1 "<<T1.block<3,3>(0,0)<<std::endl;
            Eigen::Isometry3d odom = Eigen::Isometry3d::Identity();

            odom.rotate(T1.block<3,3>(0,0));
            odom.pretranslate(T1.block<3,1>(0,3));
            //std::cout<<i<<" odom: "<<odom.matrix().cast<float>()<<std::endl;

            pcl::transformPointCloud(*ImPro_._segmented_cloud_pure, *global_cloud,odom.matrix().cast<float>());
            pcl::transformPointCloud(*ImPro_._ground_cloud, *global_ground_cloud,odom.matrix().cast<float>());
            pcl::transformPointCloud(*ImPro_._full_cloud, *global_full_cloud,odom.matrix().cast<float>());
            //std::string fileName = "/home/wjk/catkin_ws_vbdfusion/results/kittiframe.pcd";
            //pcl::io::savePCDFile(fileName, *ImPro_._segmented_cloud_pure, true);
               
           //exit(0);
            //std::cout<<"global ground cloud: "<<global_ground_cloud->size()<<std::endl;
            SM_.Mapaggregate(global_cloud, global_ground_cloud,global_full_cloud,odom.matrix().cast<float>());
            const auto& x = odom.translation()[0];
            const auto& y = odom.translation()[1];
            const auto& z = odom.translation()[2];
            PointType view;
            view.x = x;
            view.y = y;
            view.z = z;
            SM_.Saveviewpoint(view);
      }
}

void ComputeNormalVector(pcl::PointCloud<PointT>::Ptr cloud_, std::vector<int> index_,
Eigen::Vector3d viewpoint,Eigen::Vector3d& c, bool use_svd)
{
     if(use_svd){
            Eigen::Matrix<double, 3, -1> neighbors(3, index_.size());
            for (int j = 0; j < index_.size(); j++) {
            neighbors.col(j) = cloud_->points[index_[j]].getVector3fMap().cast<double>();
            }
            
            Eigen::Vector3d mean_pt = neighbors.rowwise().mean().block<3,1>(0,0);
            neighbors.colwise() -= neighbors.rowwise().mean().eval();
            Eigen::Matrix3d cov = neighbors * neighbors.transpose() / (double) index_.size();  
/*
            EigenSolver<Eigen::Matrix3d> es(cov);
            c(0) = es.eigenvectors()(2,0).real();
            c(1) = es.eigenvectors()(2,1).real();
            c(2) = es.eigenvectors()(2,2).real();
            c.normalize();
            */
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
            
            c = svd.matrixU().col(2);
            c.normalize();
            if(c.dot(viewpoint-mean_pt)<0)
                c = -1 *c;
            return;
      }


    //const int point_num = 20;
    Eigen::Matrix<double, Dynamic,3> A(index_.size(),3);
    Eigen::Matrix<double, Dynamic,1>  b(index_.size(),1);
    b.setOnes();
    b *= -1.0f;
   for(int i = 0; i < index_.size(); i++)
   {
           A(i, 0) = cloud_->points[index_[i]].x;
           A(i, 1) = cloud_->points[index_[i]].y;
           A(i, 2) = cloud_->points[index_[i]].z;
   }
    /*
    const int point_num = 20;
    Eigen::Matrix<double, point_num,3> A;
    Eigen::Matrix<double, point_num,1>  b;
    b.setOnes();
    b *= -1.0f;
   for(int i = 0; i < index_.size(); i++)
   {
           A(i, 0) = cloud_->points[index_[i]].x;
           A(i, 1) = cloud_->points[index_[i]].y;
           A(i, 2) = cloud_->points[index_[i]].z;
   }*/
   //std::cout<<"j "<<j<<std::endl;

    c = A.colPivHouseholderQr().solve(b);
    c.normalize();
    if(c.dot(viewpoint-cloud_->points[index_[0]].getVector3fMap().cast<double>())<0)
         c = -1 *c;
        
     //std::cout<<"U: "<<svd.matrixU().col(2).transpose()<<std::endl;
     //return eigenSolver.eigenvectors().col(2);
    //return svd.matrixU().col(2);
    //auto values = Eigen::Vector3d(1, 1, 1e-3);
    //Eigen::Matrix3d cov_ = svd.matrixU() * values.asDiagonal() * svd.matrixV().transpose();
}


void vdbfusion::VDBVolumeNode::checkDominateCoor(Eigen::Vector3d normvec, Eigen::Vector3d& a,Eigen::Vector3d& b,Eigen::Vector3d& c)
{
       int ind;
       normvec.normalize();

       double max_value  = normvec.cwiseAbs().minCoeff(&ind);
       //std::cout<<"normvec: "<<normvec.transpose()<<" "<<*ind<<std::endl;
        c = normvec;
        Eigen::Vector3d Cx(1,0,0);
        Eigen::Vector3d Cy(0,1,0);
        Eigen::Vector3d Cz(0,0,1);
       if(ind==0)
       {
            a = normvec.cross(Cx);
            a.normalize();
            b = normvec.cross(a);
            b.normalize();
            return;
       }
       if(ind==1)
       {
            a = normvec.cross(Cy);
            a.normalize();
            b = normvec.cross(a);
            b.normalize();
            return;
       }
       if(ind==2)
       {
            a = normvec.cross(Cz);
            a.normalize();
            b = normvec.cross(a);
            b.normalize();
            return;
       }
}

void vdbfusion::VDBVolumeNode::MapIntegrate()
{
                //  std::cout << "\033[?25l";
               pcl::PointCloud<PointT>::Ptr frame_cloud(new pcl::PointCloud<PointT>());
               pcl::PointCloud<PointT>::Ptr global_ground_cloud(new pcl::PointCloud<PointT>());
               pcl::PointCloud<PointT>::Ptr frame_cloud_transformed(new pcl::PointCloud<PointT>());
                pcl::PointCloud<PointT>::Ptr global_map_step(new pcl::PointCloud<PointT>());
                pcl::PointCloud<PointT>::Ptr global_map_step_transformed(new pcl::PointCloud<PointT>());
                pcl::PointCloud<PointT>::Ptr local_map_for_sc(new pcl::PointCloud<PointT>());
               
                int box = 3;
                float g_res = 0.05;
                std::vector<Eigen::Vector3d> point_set, origin_set;
                std::vector<float> normal_weight;
                std::vector<float> dist_weight;
                point_set.clear();
                origin_set.clear();
                normal_weight.clear();
                dist_weight.clear();
                std::vector<int> indexs;
                std::vector<float> dists;
               

                pcl::ApproximateVoxelGrid<PointType> sor,sor0,sor1,sor2;
                pcl::PointCloud<PointType>::Ptr global_ground_cloud_filtered(new pcl::PointCloud<PointT>() );
                pcl::PointCloud<PointType>::Ptr global_map_cloud_filtered(new pcl::PointCloud<PointT>() );
                 pcl::PointCloud<PointType>::Ptr ground_cloud_filtered(new pcl::PointCloud<PointT>() );
                 pcl::PointCloud<PointType>::Ptr seg_cloud_filtered(new pcl::PointCloud<PointT>() );
                  pcl::PointCloud<PointType>::Ptr full_cloud_ground_(new pcl::PointCloud<PointT>() );

              // bool test_vdb = false;
               if(test_vdb)
               {
                   for(int i = 0; i < SM_. full_cloud_set_.size(); i++){
                          auto full_cloud = SM_. full_lidar_scans_.at(i).makeShared();
                          auto scan = pcl2SensorMsgToEigen(full_cloud);
                          Eigen::Vector3d origin = SM_.viewpointset.at(i).getVector3fMap().cast<double>();
                          if(i%interval_==0 and test_ours_)
                          {
                                *Test_global_map_ += *full_cloud;
                                continue;
                           }
                          //   clock_t t3 = clock();
                          vdb_volume_.Integrate(scan, origin, [](float /*unused*/) { return 1.0; });
                   }
                   sor.setInputCloud(Test_global_map_);
                    sor.setLeafSize(0.05f, 0.05f, 0.05f); //0.1f
                        // sor.setDownsampleAllData(true);//
                    sor.filter(*Test_global_map_);
                   return;
               }
  
               point_set.clear();
               origin_set.clear();
               normal_weight.clear();
               //tree_->setInputCloud(global_map_cloud_filtered);
               
                for(int i = 0; i < SM_. full_cloud_set_.size(); i++){
                        std::cout<<"processing   i: "<<i<<std::endl;
                        auto full_cloud = SM_. full_cloud_set_.at(i).makeShared();
                        if(i>0 and i%100==0)
                        {
                            std::cout<<i/100<<" start ground integrate "<<point_set.size()<<" "<<origin_set.size()<<std::endl;
                            vdb_volume_.IntegrateVoxel(point_set, origin_set,  [](float /*unused*/) { return 1.0; }); //spave_curve_dist,
                            std::cout<<"finish ground integrate"<<std::endl;
                            point_set.clear();
                            origin_set.clear();
                        }
                        if(i % interval_==0 and test_ours_)
                        {
                            *Test_global_map_ += *SM_. full_lidar_scans_.at(i).makeShared();
                            continue;
                        }
             
                        // auto frame_cloud = SM_.lidar_scans_seg_.at(i).makeShared();
                        auto frame_cloud = SM_.full_lidar_scans_.at(i).makeShared();
                        auto scan = pcl2SensorMsgToEigen(frame_cloud);
                        Eigen::Vector3d origin = SM_.viewpointset.at(i).getVector3fMap().cast<double>();
                        //   clock_t t3 = clock();
                        vdb_volume_.Integrate(scan, origin, [](float /*unused*/) { return 1.0; });
                        //continue; // if you want to run the raw vdbfusion, uncomment this.
                        std::cout << "\033[?25l";
                        
                        for (size_t m = 0; m < ImPro_._vertical_scans-1; ++m) {
                            #pragma omp parallel for num_threads(16) schedule(guided, 8)
                            for (size_t n = 0; n < ImPro_._horizontal_scans; ++n) {
                                PointType pt = full_cloud->points[n + m * ImPro_._horizontal_scans];
                                PointType pt0 = full_cloud->points[n +2  + m * ImPro_._horizontal_scans];
                                PointType pt1 = full_cloud->points[n + (m+1) * ImPro_._horizontal_scans];
                                //std::cout<<"intensity: "<<pt.intensity<<" "<<pt0.intensity<<" "<<pt1.intensity<<std::endl;
                            if (pt.intensity >= 1 and pt0.intensity >=1  and pt1.intensity >=1)
                            {
                                std::cout<<"m n "<<m<<" "<<n<<" "<<"\r";
                                
                                 pcl::PointCloud<PointType>::Ptr agumented_cloud;
                                 if(pt.intensity >1 and pt0.intensity >1  and pt1.intensity >1)
                                     agumented_cloud = TriUpSample(pt, pt0, pt1, res_);
                                 else
                                      agumented_cloud = TriUpSample(pt, pt0, pt1, 0.5*res_);
                                 

                                Eigen::Vector3f  v0 = pt0.getVector3fMap() - pt.getVector3fMap();
                                 Eigen::Vector3f  v1 = pt1.getVector3fMap() - pt.getVector3fMap();
                                v0.normalize();
                                v1.normalize();
                                 Eigen::Vector3f v2 =  v1.cross(v0);
                                v2.normalize();

                                 Eigen::Vector3d dir = origin - pt.getVector3fMap().cast<double>();
                                 dir.normalize();
                                 if(abs(dir.dot(v2.cast<double>()))<0.15 and angle_filter_)
                                     continue;

                                // *full_cloud_ground_ += *agumented_cloud;
                                for(int j = 0; j < agumented_cloud->size(); ++j)
                                {
                                     PointType new_pt = agumented_cloud->points[j];

                                     float loc_xyz[3];
                                     loc_xyz[0] = new_pt.x / g_res;
                                     loc_xyz[1] = new_pt.y / g_res;
                                     loc_xyz[2] = new_pt.z / g_res;
                                    if (loc_xyz[0] < 0) {
                                                            loc_xyz[0] -= 1.0;
                                                        }
                                    if (loc_xyz[1] < 0) {
                                                            loc_xyz[1] -= 1.0;
                                                        }
                                    if (loc_xyz[2] < 0) {
                                                            loc_xyz[2] -= 1.0;
                                                        }

                                    VOXEL_LOC position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1],
                                                                        (int64_t)loc_xyz[2]);
                                    auto iter = voxel_map_.find(position);
                                    if (iter == voxel_map_.end())
                                    {
                                        //this is virtual viewpoint augmentation
                                        point_set.push_back(new_pt.getVector3fMap().cast<double>());
                                        origin_set.push_back(new_pt.getVector3fMap().cast<double>() + v2.cast<double>());
                                        /*
                                        for(int h  = 0; h < 3; ++h)
                                        {
                                                     point_set.push_back(new_pt.getVector3fMap().cast<double>());
                                                     origin_set.push_back(agumented_cloud->points[h].getVector3fMap().cast<double>() + v2.cast<double>());
                                        }*/
                                              //this is lidar points augmentation  
                                        point_set.push_back(new_pt.getVector3fMap().cast<double>());
                                        origin_set.push_back(origin);
                                        voxel_map_[position] = new_pt;
                                    }
                                     
                                }
                                 auto pt2 = full_cloud->points[n+2 + (m+1) * ImPro_._horizontal_scans];
                                if((int)pt2.intensity == (int)pt0.intensity)
                                {
                                        pcl::PointCloud<PointType>::Ptr agumented_cloud2;
                                        if(pt2.intensity >1 and pt0.intensity >1  and pt1.intensity >1)
                                            agumented_cloud2 = TriUpSample(pt2, pt0, pt1, res_);
                                        else
                                            agumented_cloud2 = TriUpSample(pt2, pt0, pt1, 0.5*res_);
                                        //auto agumented_cloud2 =  TriUpSample(pt2, pt0, pt1, res_);
                                        Eigen::Vector3f  vec0 = pt1.getVector3fMap() - pt2.getVector3fMap();
                                        Eigen::Vector3f  vec1 = pt0.getVector3fMap() - pt2.getVector3fMap();
                                        vec0.normalize();
                                        vec1.normalize();
                                        Eigen::Vector3f  vec2 =  vec1.cross(vec0);
                                        vec2.normalize();
                                         
                                         dir = origin - pt2.getVector3fMap().cast<double>();
                                        dir.normalize();
                                        if(abs(dir.dot(vec2.cast<double>()))<0.15 and angle_filter_)
                                            continue;

                                        //*full_cloud_ground_ += *agumented_cloud2;
                                        for(int j = 0; j < agumented_cloud2->size(); ++j)
                                        {
                                            PointType new_pt = agumented_cloud2->points[j];
                                            float loc_xyz[3];
                                            loc_xyz[0] = new_pt.x / g_res;
                                            loc_xyz[1] = new_pt.y / g_res;
                                            loc_xyz[2] = new_pt.z / g_res;
                                            if (loc_xyz[0] < 0) {
                                                                    loc_xyz[0] -= 1.0;
                                                                }
                                            if (loc_xyz[1] < 0) {
                                                                    loc_xyz[1] -= 1.0;
                                                                }
                                            if (loc_xyz[2] < 0) {
                                                                    loc_xyz[2] -= 1.0;
                                                                }

                                            VOXEL_LOC position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1],
                                                                                (int64_t)loc_xyz[2]);
                                            auto iter = voxel_map_.find(position);
                                             if (iter == voxel_map_.end())
                                            {
                                                /*
                                                for(int h  = 0; h <3; ++h)
                                                {
                                                     point_set.push_back(new_pt.getVector3fMap().cast<double>());
                                                     origin_set.push_back(agumented_cloud2->points[h].getVector3fMap().cast<double>() + vec2.cast<double>());
                                                }*/

                                                point_set.push_back(new_pt.getVector3fMap().cast<double>());
                                                origin_set.push_back(new_pt.getVector3fMap().cast<double>() + vec2.cast<double>());
                                                point_set.push_back(new_pt.getVector3fMap().cast<double>());
                                                origin_set.push_back(origin);
                                                voxel_map_[position] = new_pt;
                                            }                                            
                                        }
                                }
                                

                            }
                            }
                        }

                        

               }
                vdb_volume_.IntegrateVoxel(point_set, origin_set,  [](float /*unused*/) { return 1.0; }); //spave_curve_dist,
                std::cout<<"Test_global_map "<<Test_global_map_->size()<<std::endl;
                sor.setInputCloud(Test_global_map_);
                sor.setLeafSize(0.05f, 0.05f, 0.05f); //0.1f
                    // sor.setDownsampleAllData(true);//
                sor.filter(*Test_global_map_);

}

void vdbfusion::VDBVolumeNode::testbaselines()
{
     pcl::PolygonMesh mesh;
     //
     
     //pcl::io::loadPLYFile("/home/wjk/SLAMesh_raw_ws/kitti07voxfiled.ply",mesh);
     pcl::io::loadPLYFile("/home/wjk/Dataset/mesh_new1.ply",mesh); //please use your own data path
     //pcl::io::loadPLYFile("/home/wjk/SLAMesh_raw_ws/kitti07_testslamesh_raw_mesh.ply",mesh);
     std::cout << "Polygon: " << mesh.polygons.size() << std::endl;
    std::cout <<"cloud: "<< mesh.cloud.width << std::endl;
    
    pcl::PointCloud<PointType>::Ptr   cloud(new pcl::PointCloud<PointType>);
    fromPCLPointCloud2(mesh.cloud, *cloud);
     
    pcl::PointCloud<PointT>::Ptr global_map(new pcl::PointCloud<PointT>());
    double S = 0;
    std::vector<Eigen::Vector3d> normals_set;
    int num = 0;
     for(size_t i = 0; i < mesh.polygons.size(); i++)
     {
         
          PointT A = cloud->points.at(mesh.polygons[i].vertices[0]);
          PointT B = cloud->points.at(mesh.polygons[i].vertices[1]);
          PointT C = cloud->points.at(mesh.polygons[i].vertices[2]);
         
         PointT D;
         D.x = (A.x+B.x+C.x)/3;
         D.y = (A.y+B.y+C.y)/3;
         D.z = (A.z+B.z+C.z)/3;
         
         float a = (A.getVector3fMap() - B.getVector3fMap()).norm();
         float b = (A.getVector3fMap() - C.getVector3fMap()).norm();
         float c = (C.getVector3fMap() - B.getVector3fMap()).norm();
          
         float p = 0.5*(a+b+c);

         D.intensity = a*b*c/(4*sqrt(p*(p-a)*(p-b)*(p-c)));

          if(isnan(D.getVector3fMap().norm()) or D.getVector3fMap().norm()>1000000)
              continue;

          global_map->points.push_back(D);


          Eigen::Vector3d an =  (A.getVector3fMap() - B.getVector3fMap()).cast<double>();
          Eigen::Vector3d bn = (A.getVector3fMap() - C.getVector3fMap()).cast<double>();
          Eigen::Vector3d n = an.cross(bn);
          S += n.norm();
         num = num+1;
          normals_set.push_back(n.normalized());
     }
     global_map->width = 1;
     global_map->height  = num;
     float score = 0;//; =  computeMapdistance(global_map0);
    tree_->setInputCloud(global_map);
     std::vector<int> nn_indices;
    std::vector<float> nn_dists;
    float fitness_score = 0;
    int nr = 0;
    for(size_t i = 0; i < Test_global_map_->size(); ++i) {
         PointType pt = Test_global_map_->points[i];
         if(isnan(pt.x) or isnan(pt.y) or isnan(pt.z))
            continue;
         tree_->nearestKSearch(pt, 3, nn_indices, nn_dists);
         auto v1 = pt.getVector3fMap() - global_map->points[nn_indices[0]].getVector3fMap();
         auto n1 = normals_set.at(nn_indices[0]);
         float d1 = v1.cast<double>().dot(n1);
         float d2 = sqrt(pow(v1.norm(),2)-pow(d1,2));
         if(d2<global_map->points[nn_indices[0]].intensity and nn_dists[0]<0.5)
         {
             fitness_score += pow(abs(d1),1);
            nr++;
         }


    }

     std::cout<<"the score of baselines is "<<(float)fitness_score/nr<<" "<<S<<" "<<(float)nr/Test_global_map_->size()<<" "<<nr<<std::endl;

}

pcl::PointCloud<PointT>::Ptr vdbfusion::VDBVolumeNode::SampleGlobalMap(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    pcl::PointCloud<PointT>::Ptr global_map(new pcl::PointCloud<PointT>());
    double S = 0;
    std::vector<Eigen::Vector3d> normals_set;
    int num = 0;
     for(size_t i = 0; i < F.rows(); i++)
     {
          PointT A = point2pcl(V.row(F.row(i)[0]));
          PointT B = point2pcl(V.row(F.row(i)[1]));
          PointT C = point2pcl(V.row(F.row(i)[2]));
         
         PointT D;
         D.x = (A.x+B.x+C.x)/3;
         D.y = (A.y+B.y+C.y)/3;
         D.z = (A.z+B.z+C.z)/3;
         
         float a = (A.getVector3fMap() - B.getVector3fMap()).norm();
         float b = (A.getVector3fMap() - C.getVector3fMap()).norm();
         float c = (C.getVector3fMap() - B.getVector3fMap()).norm();
          
         float p = 0.5*(a+b+c);

         D.intensity = a*b*c/(4*sqrt(p*(p-a)*(p-b)*(p-c)));

          if(isnan(D.getVector3fMap().norm()) or D.getVector3fMap().norm()>1000000)
              continue;

          global_map->points.push_back(D);


          Eigen::Vector3d an =  (A.getVector3fMap() - B.getVector3fMap()).cast<double>();
          Eigen::Vector3d bn = (A.getVector3fMap() - C.getVector3fMap()).cast<double>();
          Eigen::Vector3d n = an.cross(bn);
          S += n.norm();
         num = num+1;
          normals_set.push_back(n.normalized());
     }
     global_map->width = 1;
     global_map->height  = num;
     float score = 0;//; =  computeMapdistance(global_map0);
    tree_->setInputCloud(global_map);
     std::vector<int> nn_indices;
    std::vector<float> nn_dists;
    float fitness_score = 0;
    int nr = 0;
    for(size_t i = 0; i < Test_global_map_->size(); ++i) {
         PointType pt = Test_global_map_->points[i];
         if(isnan(pt.x) or isnan(pt.y) or isnan(pt.z))
            continue;
         tree_->nearestKSearch(pt, 3, nn_indices, nn_dists);
         auto v1 = pt.getVector3fMap() - global_map->points[nn_indices[0]].getVector3fMap();
         auto n1 = normals_set.at(nn_indices[0]);
         float d1 = v1.cast<double>().dot(n1);
         float d2 = sqrt(pow(v1.norm(),2)-pow(d1,2));
         if(d2<global_map->points[nn_indices[0]].intensity and nn_dists[0]<0.5)
         {
             fitness_score += pow(abs(d1),1);
            nr++;
         }


    }

     std::cout<<"the score is "<<(float)fitness_score/nr<<" "<<S<<" "<<(float)nr/Test_global_map_->size()<<" "<<nr<<std::endl;
     std::string n = "/home/wjk/catkin_ws_vbdfusion/results/globalpcdfromply.pcd";
     pcl::io::savePCDFileASCII(n, *global_map); //将点云保存到PCD文件中
     std::string n0 = "/home/wjk/catkin_ws_vbdfusion/results/testglobalmap.pcd";
     pcl::io::savePCDFileASCII(n0, *Test_global_map_); //将点云保存到PCD文件中
     return global_map;
}

bool vdbfusion::VDBVolumeNode::saveVDBVolume(vdbfusion_ros::save_vdb_volume::Request& path,
                                             vdbfusion_ros::save_vdb_volume::Response& response) {
        ROS_INFO("Reading Pose from pose.txt");
         
        kittipose_path = "/home/wjk/Dataset/kitti/dataset/poses/07.txt";
        //kittipose_path = "/home/wjk/Dataset/MaiCity/mai_city/bin/poses/01.txt";
        load_pose_kitti();
        
        std::cout<<transpose.at(0)<<std::endl;
        //std::cout<<"pose size: "<<transpose.size()<<std::endl;
        //exit(0);
        Integrate();
         
       
        MapIntegrate();


        ROS_INFO("Generating the mesh and VDB grid files ...");
        
    ROS_INFO("Saving the mesh and VDB grid files ...");
    //SM_.Mapvisualize();
    std::cout<<"num of scans: "<<cloud_integrate_num_<<std::endl;
    std::string volume_name = path.path;
    openvdb::io::File(volume_name + "_grid.vdb").write({vdb_volume_.tsdf_});

    // Run marching cubes and save a .ply file
    auto [vertices, triangles] =
        this->vdb_volume_.ExtractTriangleMesh(this->fill_holes_, this->min_weight_);

    Eigen::MatrixXd V(vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); i++) {
        V.row(i) = Eigen::VectorXd::Map(&vertices[i][0], vertices[i].size());
    }

    Eigen::MatrixXi F(triangles.size(), 3);
    for (size_t i = 0; i < triangles.size(); i++) {
        F.row(i) = Eigen::VectorXi::Map(&triangles[i][0], triangles[i].size());
    }
    std::cout<<"F.row(0) "<<F.row(0)<<std::endl<<F.row(1)<<std::endl;
    std::cout<<"V.row(0) "<<V.row(0)<<std::endl<<V.row(1)<<std::endl;
    igl::write_triangle_mesh(volume_name + "_mesh.ply", V, F, igl::FileEncoding::Binary);
    //SampleGlobalMap(V, F);
    //ROS_INFO("Computing the Chamfer Distance");
    //ComputeChamferDistance(V);
    ROS_INFO("Computing the TSDF Error");
    //ComputeTSDFError(); //No need for kitti

    SampleGlobalMap(V, F);
    ROS_INFO("Done saving the mesh and VDB grid files");
    testbaselines();
    //ROS_INFO("Done testing baselines");

    
    return true;
}

int main(int argc, char** argv) {
    ros::init(argc, argv, "vdbfusion_rosnode");
    vdbfusion::VDBVolumeNode vdb_volume_node;
    ros::spin();
}
