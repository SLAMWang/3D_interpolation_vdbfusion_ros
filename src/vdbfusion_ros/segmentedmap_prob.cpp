// Copyright 2013, Ji Zhang, Carnegie Mellon University
// Further contributions copyright (c) 2016, Southwest Research Institute
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "segmentedmap.h"
#define NUM_MATCH_POINTS (3)

bool SegmentedMap::esti_plane(Eigen::Vector4f &pca_result, std::vector<PointType> points) {
  Eigen::MatrixXf A(points.size(),3);
  Eigen::MatrixXf  b(points.size(),1);
  b.setOnes();
  b *= -1.0f;

  for (int j = 0; j < points.size(); j++) {
    A(j, 0) = points[j].x;
    A(j, 1) = points[j].y;
    A(j, 2) = points[j].z;
  }

  Eigen::Matrix<float, 3, 1> normvec = A.colPivHouseholderQr().solve(b);



  pca_result(0) = normvec(0) ;
  pca_result(1) = normvec(1);
  pca_result(2) = normvec(2);
  pca_result(3) = 1.0;
  float mean_dist = 0;
   for(int i = 0; i < points.size(); ++i)
   {
        PointType pt = points.at(i);
        mean_dist = mean_dist + fabs(pt.x *pca_result(0)+pt.y *pca_result(1)+pt.z *pca_result(2)+1);
   }
   mean_dist = mean_dist/(float)points.size();
  //std::cout<<"esti_plane"<<std::endl;
   if(mean_dist > 0.1)
      return false;
   else
      return true;

}


float SegmentedMap::computeSegdistance(pcl::PointCloud<PointType>::Ptr cloudtarget,pcl::PointCloud<PointType>::Ptr cloudsource)
{
   
     /*
     Eigen::Matrix<double, 4, -1> neighbors(4, cloudtarget->size());
    for (int j = 0; j < cloudtarget->size(); j++) {
      neighbors.col(j) = cloudtarget->at(j).getVector4fMap().template cast<double>();
    }
    
    Eigen::Vector3d mean_pt = neighbors.rowwise().mean().block<3,1>(0,0);
    neighbors.colwise() -= neighbors.rowwise().mean().eval();
    Eigen::Matrix4d cov = neighbors * neighbors.transpose() / (double)cloudtarget->size();  
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(cov.block<3, 3>(0, 0), Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    auto values = Eigen::Vector3d(1, 1, 1e-3);
    Eigen::Matrix3d cov_ = svd.matrixU() * values.asDiagonal() * svd.matrixV().transpose();
    
    double fitnessscore = 0;
    int nr  =  0;
     for(size_t i = 0; i < cloudsource->size(); ++i) {
         PointType pt = cloudsource->points[i];
         Eigen::Vector3d point = pt.getVector3fMap().template cast<double>();
         Eigen::Vector3d dist = point - mean_pt;
         double d = std::sqrt(dist.transpose() * cov_ *dist);
         if(d < 0.25)
         {
             fitnessscore  += d;
             nr++;
         }
     }
     if(nr>0.1*cloudsource->size())
        return (float)fitnessscore/nr;
     else
         return 100;
*/
   

  /*
    pcl::search::KdTree<PointType>::Ptr tree_(new pcl::search::KdTree<PointType>());
      tree_->setInputCloud(cloudtarget);

      double fitness_score = 0.0;

      std::vector<int> nn_indices(1);
      std::vector<float> nn_dists(1);
      // For each point in the source dataset
      int nr = 0;
      int point_num = 3;
      for(size_t i = 0; i < cloudsource->size(); ++i) {
        // Find its nearest neighbor in the target
        tree_->nearestKSearch(cloudsource->points[i], point_num, nn_indices, nn_dists);
          // Add to the fitness score
        
        if(nn_dists[0]>0.5)
             continue;

         std::vector<PointType> points;
         for(int k = 0; k < point_num; k++)
              points.push_back(cloudtarget->points[nn_indices[k]]);
        Eigen::Vector4f pca_result;
        esti_plane(pca_result, points);
        fitness_score =  fabs(pca_result(0) * cloudsource->points[i].x + pca_result(1) * cloudsource->points[i].y +
             pca_result(2) * cloudsource->points[i].z + pca_result(3));
        
        if(fitness_score<0.05)
              nr++;
      }
      float r = (float)nr/cloudsource->size();
      return 1-r;*/
  
  pcl::search::KdTree<PointType>::Ptr tree_(new pcl::search::KdTree<PointType>());
  tree_->setInputCloud(cloudtarget);

  double fitness_score = 0.0;

  std::vector<int> nn_indices(1);
  std::vector<float> nn_dists(1);

  
  // For each point in the source dataset
  int nr = 0;
  int point_num = 1;
  
  for(size_t i = 0; i < cloudsource->size(); ++i) {
    // Find its nearest neighbor in the target
    PointType pt = cloudsource->points[i];
    tree_->nearestKSearch(pt, point_num, nn_indices, nn_dists);
    
      // Add to the fitness score
       if( nn_dists[0] < 0.25)
         {
             fitness_score  += nn_dists[0]/std::sqrt(pow(pt.x,2)+pow(pt.y,2)+pow(pt.z,2));
             nr++;
         }

  }
     if(nr>0.1*cloudsource->size())
        return (float)fitness_score/nr;
     else
         return 100;
}


void computeSVD(Eigen::Vector3f& center, Eigen::Matrix3f&  cov, float& pdist, pcl::PointCloud<PointType>::Ptr s_cloud)
{
    Eigen::Matrix<float, 3, -1> neighbors(3, s_cloud->size());
    for (int j = 0; j < s_cloud->size(); j++) {
      neighbors.col(j) = s_cloud->at(j).getVector3fMap();
    }
    
    center = neighbors.rowwise().mean();
    neighbors.colwise() -= neighbors.rowwise().mean().eval();
    cov = neighbors * neighbors.transpose() / (float)s_cloud->size();
    //Eigen::JacobiSVD<Eigen::Matrix3f> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //Eigen::Vector3f values = Eigen::Vector3f(1, 1, 1e-3);
      
    //cov = svd.matrixU() * values.asDiagonal() * svd.matrixV().transpose();
     pdist = 0;
 
}

bool SegmentedMap::SegmentedFrame(pcl::PointCloud<PointType>::Ptr cloud)
{
      //1. get the size of the segments
      int num_seg = 0;
      for(int i = 0; i < cloud->size(); ++i)
        {
            PointType pt = cloud->points[i];
            int seg_index = pt.intensity;
            if(seg_index>num_seg)
               num_seg = seg_index;
        }
        //curr_segment_clouds.resize(num_seg+1);
        segmented_frame_.clear();
        for(int i = 0; i < num_seg+1; i++)
        {
            pcl::PointCloud<PointType>::Ptr pc;
            pc.reset(new pcl::PointCloud<PointType>());
             segmented_frame_.push_back(pc);
        }
        int max_index = 0;
        std::vector<int> cloud_size;
        cloud_size.resize(num_seg+1,0);
        for(int i = 0; i < cloud->size(); ++i)
        {
            PointType pt = cloud->points[i];
            int seg_index = pt.intensity;
            pt.intensity = frame_index_;
            segmented_frame_.at(seg_index)->points.push_back(pt);
            int a = cloud_size.at(seg_index);
            a++;
            cloud_size.at(seg_index)=a;
        }
        //std::cout<<"frame segment number "<<num_seg<<std::endl;

        for(int i = 0; i < cloud_size.size();++i)
        {
            //std::cout<<"size: "<<segmented_frame_.at(i)->size()<<" "<<cloud_size.at(i)<<std::endl;
            segmented_frame_.at(i)->width = 1;
            segmented_frame_.at(i)->height = cloud_size.at(i);
        }

        Mapmerge();
        return true;
}

void SegmentedMap::Mapmerge(){
  //Then, merge to the global map
        //std::cout<<"segment_clouds_ "<<segmented_map_.size()<<std::endl;
        if(segmented_map_.size()==0)
        {
            //it is the first scan 

            for(int i = 0; i < segmented_frame_.size(); ++i)
            {
              if(segmented_frame_.at(i)->size()<10)
                  continue;
                Seg_ segment_;
                segment_.seg_cloud_.reset(new pcl::PointCloud<PointType>());
                pcl::copyPointCloud(*segmented_frame_.at(i), *segment_.seg_cloud_);
                Eigen::Vector3f center;
                Eigen::Matrix3f cov;
                float pdist;
                computeSVD(center,cov,pdist,segment_.seg_cloud_);
                if(std::max(cov(0,0),std::max(cov(1,1),cov(2,2)))<0.01)
               {
                   //continue;
                   //std::string n = "/home/wjk/catkin_ws_vbdfusion/results/smallsegment.pcd";
                   //pcl::io::savePCDFileASCII(n, *segment_.seg_cloud_); //将点云保存到PCD文件中
               }
                segment_.merged_center_ = center;
                segment_.merged_cov_ = cov;
                segment_.last_center_ = center;
                segment_.last_cov_ = cov;
                //std::cout<<"center: "<<center(0)<<" "<<center(1)<<" "<<center(2)<<" "<<std::endl;
                //std::cout<<"cov: "<<cov(0,0)<<" "<<cov(1,1)<<" "<<cov(2,2)<<" "<<segment_.seg_cloud_->size()<<std::endl;
                segment_.obs_times_ = 1; 
                segment_.dist_ = pdist;
                segmented_map_.push_back(segment_);
               // std::cout<<"pdist: "<<pdist<<" "<<segment_.seg_cloud_->size()<<std::endl;
               
              

            }
            //std::cout<<"init map"<<std::endl;
            return;
        }
        for(int i = 0; i < segmented_frame_.size(); ++i)
        {
             if(segmented_frame_.at(i)->size()<10)
                  continue;

             Eigen::Vector3f center;
             Eigen::Matrix3f cov;
             float pdist;
            computeSVD(center,cov,pdist,segmented_frame_.at(i));
            //std::cout<<"seg_clouds_ "<<seg_cloud.size()<<" "<<segmented_frame_.at(i)->size()<<std::endl;
             if(std::max(cov(0,0),std::max(cov(1,1),cov(2,2)))<0.01)
               {
                continue;
                   //std::string n = "/home/wjk/catkin_ws_vbdfusion/results/smallsegment.pcd";
                   //pcl::io::savePCDFileASCII(n, *segment_.seg_cloud_); //将点云保存到PCD文件中
               }
            float min_score = 100000000;
            float score;
            int min_index;
            float min_last_dist;
            for(int j = 0; j < segmented_map_.size();++j)
            {
                Seg_ curr_mapseg = segmented_map_.at(j);
                
                Eigen::Vector3f ptn = center-curr_mapseg.last_center_;
                if(ptn.norm()>3)
                     continue;
                
                //std::cout<<"ptn: "<<j<<" "<<ptn.transpose()<<std::endl;
                //std::cout<<"cov: "<<cov(0,0)<<" "<<cov(1,1)<<" "<<cov(2,2)<<" "<<std::endl;
                //std::cout<<"last_cov: "<<curr_mapseg.last_cov_(0,0)<<" "<<curr_mapseg.last_cov_(1,1)<<" "<<curr_mapseg.last_cov_(2,2)<<" "<<std::endl;

                score= ptn.transpose()*(cov+curr_mapseg.last_cov_).inverse()*ptn;//
                //*ptn;
                //score  = exp(-1*pow(ptn(0),2)/(curr_mapseg.last_cov_(0,0)+cov(0,0)))+
                //exp(-1*pow(ptn(1),2)/(curr_mapseg.last_cov_(1,1)+cov(1,1)));
                //exp(-1*pow(ptn(2),2)/(curr_mapseg.last_cov_(2,2)+cov(2,2)));
               
                if(score<min_score)
                {
                    min_score = score;
                    min_index = j;
                    min_last_dist = curr_mapseg.dist_;
                }
            }
            //std::cout<<"min_score "<<min_score<<" "<<min_index<<" "<<min_last_dist<<std::endl;
             if(min_score<1)
            {
                  *segmented_map_.at(min_index).seg_cloud_ += *segmented_frame_.at(i);
                  segmented_map_.at(min_index).last_center_ = center;
                  segmented_map_.at(min_index).last_cov_ = cov;
                  //computeSVD(segmented_map_.at(min_index).merged_center_,segmented_map_.at(min_index).merged_cov_,
                  //segmented_map_.at(min_index).dist_,segmented_map_.at(min_index).seg_cloud_);

                  segmented_map_.at(min_index).obs_times_++;
            }else{
                   Seg_ segment_;
                   segment_.merged_center_ = center;
                   segment_.merged_cov_ = cov;
                    segment_.last_center_ = center;
                    segment_.last_cov_ = cov;
                    segment_.obs_times_ = 1; 
                    segment_.seg_cloud_.reset(new pcl::PointCloud<PointType>());
                    pcl::copyPointCloud(*segmented_frame_.at(i), *segment_.seg_cloud_);
                  segmented_map_.push_back(segment_);
            }              
        }
}

void SegmentedMap::Mapweighted()
{
      //std::cout<<"weighted"<<std::endl;
      kdtree_->setInputCloud(global_map_);
      float dist = 0.5;
      std::vector<int> nn_indices(1);
      std::vector<float> nn_dists(1);
      Eigen::Vector4f pca_result;
      std::vector<PointType> neighbor_points;
      float pt_2_plane;
      for(int i = 0; i < current_weighted_cloud_->size(); ++i)
      {
          PointType pt = current_weighted_cloud_->points[i];
           kdtree_->nearestKSearch(pt, 3, nn_indices, nn_dists);
           neighbor_points.clear();
           
           for(int j = 0; j < nn_indices.size(); ++j)
                 neighbor_points.push_back(global_map_->points[nn_indices[j]]);
          
          if(esti_plane(pca_result,neighbor_points))
          {
              if(fabs(pca_result(0)*pt.x + pca_result(1)*pt.y + pca_result(2)*pt.z+1)<0.2)
              {
                 for(int j = 0; j < nn_indices.size(); ++j)
                     global_map_->points[nn_indices[j]].intensity++;
              }
          }
        }
}

void SegmentedMap::Mapaggregate(pcl::PointCloud<PointType>::Ptr cloud, Eigen::Matrix4f odom){

  //pcl::copyPointCloud(*cloud,*current_weighted_cloud_);
  pcl::PointCloud<PointType> seg_cloud = *cloud;
  //SegmentedFrame(seg_cloud.makeShared());
  lidar_scans_seg_.push_back(seg_cloud);
   for(int i = 0; i < cloud->size(); i++)
   {
       cloud->points[i].intensity = frame_index_;
   }
    *global_map_ += *cloud;
    frame_index_++;
    lidar_scans_.push_back(*cloud);
    scan_odoms_.push_back(odom);
    //Mapweighted();
    //SegmentedFrame(cloud);
    
}


void SegmentedMap::Saveviewpoint(PointType p)
{
   viewpointset.push_back(p);
}

void SegmentedMap::Mapvisualize()
{
  /*
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr rgb_pointcloud(
        new pcl::PointCloud<pcl::PointXYZRGB>);
        int r0 = 250;
        int g0 = 250;
        for (std::size_t i = 0; i < global_map_->size(); ++i) {
          //generate color
               
                int b0 =  (int)global_map_->points[i].intensity%255 + 1;
                pcl::PointXYZRGB rgb_point;
                const auto &point = global_map_->points[i];
                rgb_point.x = point.x;
                rgb_point.y = point.y;
                rgb_point.z = point.z;
                rgb_point.r = r0;
                rgb_point.g = g0;
                rgb_point.b = b0;
                rgb_pointcloud->push_back(std::move(rgb_point));
        }

    std::string n = "/home/wjk/catkin_ws_vbdfusion/results/colormap.pcd";
     pcl::io::savePCDFileASCII(n, *rgb_pointcloud); //将点云保存到PCD文件中
*/
  
      pcl::PointCloud<pcl::PointXYZRGB>::Ptr rgb_pointcloud(
      new pcl::PointCloud<pcl::PointXYZRGB>);
    for (std::size_t i = 0; i < segmented_map_.size(); ++i) {
      if(segmented_map_.at(i).obs_times_>5)
          continue;
      auto seg_cloud = segmented_map_.at(i).seg_cloud_;

      *global_map_pure_ += *seg_cloud;

      //generate color
      int r0 = rand()%255 + 1;
      int g0 = rand()%255 + 1;
      int b0 = rand()%255 + 1;
       for(int j = 0; j < seg_cloud->size(); j++)
       {
            pcl::PointXYZRGB rgb_point;
            const auto &point = seg_cloud->points.at(j);
            rgb_point.x = point.x;
            rgb_point.y = point.y;
            rgb_point.z = point.z;
            rgb_point.r = r0;
            rgb_point.g = g0;
            rgb_point.b = b0;
            rgb_pointcloud->push_back(std::move(rgb_point));
       }
    }

    std::string n = "/home/wjk/catkin_ws_vbdfusion/results/colormap.pcd";
     pcl::io::savePCDFileASCII(n, *rgb_pointcloud); //将点云保存到PCD文件中

/*
    pcl::visualization::PCLVisualizer::Ptr visualizer(
        new pcl::visualization::PCLVisualizer("PointCloud Visualizer"));
    visualizer->setBackgroundColor(0, 0, 0);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>
        rgb_color_handler(rgb_pointcloud);
    visualizer->addPointCloud<pcl::PointXYZRGB>(rgb_pointcloud, rgb_color_handler,
                                                "RGB PointCloud");
    visualizer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "RGB PointCloud");
    visualizer->addCoordinateSystem(5.0);*/
}
