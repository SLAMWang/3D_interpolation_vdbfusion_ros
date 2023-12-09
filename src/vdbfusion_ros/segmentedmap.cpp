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


void computeSVD(Eigen::Vector3f& center, Eigen::Matrix3f&  cov, float& pdist,
 pcl::PointCloud<PointType>::Ptr s_cloud, Vertexes& _vert)
{
    Eigen::Matrix<float, 3, -1> neighbors(3, s_cloud->size());
    for (int j = 0; j < s_cloud->size(); j++) {
      neighbors.col(j) = s_cloud->at(j).getVector3fMap();
      //neighbors.col(j)(2) = 0;
    }
    
    center = neighbors.rowwise().mean();
    neighbors.colwise() -= neighbors.rowwise().mean().eval();
    cov = neighbors * neighbors.transpose() / (float)s_cloud->size();
    Eigen::JacobiSVD<Eigen::Matrix3f> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //Eigen::Vector3f values = Eigen::Vector3f(1, 1, 1e-3);
    
   // std::cout<<"svd.matrixU() "<<std::endl<<svd.matrixU()<<std::endl;
   // std::cout<<"eigenvalues: "<<svd.singularValues().transpose()<<std::endl;
    //cov = svd.matrixU() * values.asDiagonal() * svd.matrixV().transpose();

    float max_x = -1;
    float max_y = -1;
    float min_x = 100000;
    float min_y = 100000;
    Eigen::Vector2d eigen_v_x(svd.matrixU()(0,0),svd.matrixU()(1,0));
    Eigen::Vector2d eigen_v_y(svd.matrixU()(0,1),svd.matrixU()(1,1));
    eigen_v_x.normalized();
    eigen_v_y.normalized();
    for(int i =0; i < s_cloud->size(); i++)
    {
         PointType pt = s_cloud->points[i];
         Eigen::Vector2d pointV(pt.x-center[0],pt.y-center[1]);
         float dx = pointV.dot(eigen_v_x);
         float dy = pointV.dot(eigen_v_y);

         //std::cout<<i<< "i dx dy "<<dx<<" "<<dy<<std::endl;

         if(dx>max_x)
             max_x = dx;
         if(dy>max_y)
             max_y = dy;
          if(dx<min_x)
             min_x = dx;
          if(dy<min_y)
             min_y = dy; 
    }
    max_x += 0.05;
    max_y += 0.05;

   // std::cout<<"max_x "<<max_x<<" "<<max_y<<" "<<min_x<<" "<<min_y<<std::endl;
    Eigen::Vector2d center2d(center[0],center[1]);
    Eigen::Vector2d P0 = center2d + max_x * eigen_v_x + max_y * eigen_v_y;
    Eigen::Vector2d P1 = center2d + max_x * eigen_v_x + min_y * eigen_v_y;
    Eigen::Vector2d P2 = center2d + min_x * eigen_v_x + min_y * eigen_v_y;
    Eigen::Vector2d P3 = center2d + min_x * eigen_v_x + max_y * eigen_v_y;
    
        //_vert.resize(4);

        _vert.push_back(IOU::Point(P0[0],P0[1]));
        _vert.push_back(IOU::Point(P1[0],P1[1]));
        _vert.push_back(IOU::Point(P2[0],P2[1]));
        _vert.push_back(IOU::Point(P3[0],P3[1]));
    //beInSomeWiseEx(_vert,ClockWise);
    
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
            pt.intensity = frame_index_-1;
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
                 Vertexes vert;
                 vert.clear();
                computeSVD(center,cov,pdist,segment_.seg_cloud_,  vert);
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
                segment_._vert = vert;
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
             Vertexes vert;
             vert.clear();
            computeSVD(center,cov,pdist,segmented_frame_.at(i),vert);
            if(areaEx(vert)<0.1)
                continue;
            //std::cout<<"seg_clouds_ "<<seg_cloud.size()<<" "<<segmented_frame_.at(i)->size()<<std::endl;
             
            float min_score = 100000;
            float score;
            int min_index;
            float min_last_dist;
            for(int j = 0; j < segmented_map_.size();++j)
            {
                Seg_ curr_mapseg = segmented_map_.at(j);
                
                Eigen::Vector3f ptn = center-curr_mapseg.last_center_;
                if(ptn.norm()>5)
                     continue;
                
                //std::cout<<"ptn: "<<j<<" "<<ptn.transpose()<<std::endl;
                //std::cout<<"cov: "<<cov(0,0)<<" "<<cov(1,1)<<" "<<cov(2,2)<<" "<<std::endl;
                //std::cout<<"last_cov: "<<curr_mapseg.last_cov_(0,0)<<" "<<curr_mapseg.last_cov_(1,1)<<" "<<curr_mapseg.last_cov_(2,2)<<" "<<std::endl;
               //std::cout<<"Area: "<<vert.size()<<" "<<curr_mapseg._vert.size()<<" "<<whichWiseEx(vert)<<" "<<whichWiseEx(curr_mapseg._vert)<<" "
               //<<NoneWise<<" "<<areaEx(vert)<<" "<<areaEx(curr_mapseg._vert)<<std::endl;
               //std::cout<<areaIntersectionEx(vert,curr_mapseg._vert)<<
               //" "<<areaUnionEx(vert,curr_mapseg._vert)<<std::endl;
              
              

              if(false)
              {
                 for(int m = 0; m < curr_mapseg._vert.size(); ++m)
                 {
                    std::cout<<curr_mapseg._vert.at(m)[0]<<" "<<curr_mapseg._vert.at(m)[1]<<std::endl;
                 }
              }

              if(false)
              {
                std::cout<<"vert: "<<std::endl;
                 for(int m = 0; m < vert.size(); ++m)
                 {
                    std::cout<<vert.at(m)[0]<<" "<<vert.at(m)[1]<<std::endl;
                 }
              }

               //score = areaIntersectionEx(vert,curr_mapseg._vert)/std::min(areaEx(vert),areaEx(curr_mapseg._vert));
               //score = iouEx(vert,curr_mapseg._vert);
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
             if(min_score<20) 
            {
                  *segmented_map_.at(min_index).seg_cloud_ += *segmented_frame_.at(i);
                   vert.clear();
                 //computeSVD(center,cov,pdist,segmented_map_.at(min_index).seg_cloud_,vert);
                  segmented_map_.at(min_index).last_center_ = center;
                  segmented_map_.at(min_index).last_cov_ = cov;
                  segmented_map_.at(min_index)._vert.swap(vert);
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
                    segment_._vert.swap(vert);
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

void SegmentedMap::Mapaggregate(pcl::PointCloud<PointType>::Ptr cloud, pcl::PointCloud<PointType>::Ptr ground_cloud,
pcl::PointCloud<PointType>::Ptr full_cloud,Eigen::Matrix4f odom){
    

  //pcl::copyPointCloud(*cloud,*current_weighted_cloud_);
  pcl::PointCloud<PointType> seg_cloud = *cloud;
  //SegmentedFrame(seg_cloud.makeShared());
  lidar_scans_seg_.push_back(seg_cloud);
  full_lidar_scans_.push_back(*cloud+*ground_cloud);
  full_cloud_set_.push_back(*full_cloud);
   for(int i = 0; i < cloud->size(); i++)
   {
       cloud->points[i].intensity = frame_index_;
       
   }
   for(int i = 0; i < ground_cloud->size(); i++)
   {
       //if(isnan(ground_cloud->points[i].x) or isnan(ground_cloud->points[i].y) or isnan(ground_cloud->points[i].z))
       ground_cloud->points[i].intensity = frame_index_;
} 
/*
    pcl::NormalEstimation<PointType, pcl::Normal> ne;
	ne.setInputCloud (ground_cloud);
	pcl::search::KdTree<PointType>::Ptr tree (new pcl::search::KdTree<PointType> ());
	ne.setSearchMethod (tree);
	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
	ne.setRadiusSearch (0.5);
	//计算特征值
	ne.compute (*cloud_normals);

    for(int i = 0; i < ground_cloud->size();++i)
    {
        Eigen::Vector3d fake_normal = odom.block<3,1>(0,0).cast<double>() - ground_cloud->points.at(i).getVector3fMap().cast<double>();
        
        auto normal_vec = cloud_normals->at(i).getNormalVector3fMap().cast<double>();
        
                if(normal_vec.dot(fake_normal)<0)
                     global_ground_cloud_norm_.push_back(-1 * normal_vec);
                else
                     global_ground_cloud_norm_.push_back(normal_vec);
    }*/
  //*global_ground_cloud_normals_ += *cloud_normals;
   lidar_ground_scans_.push_back(*ground_cloud);
   *global_ground_map_ += *ground_cloud;

   
    *global_map_ += *cloud;
    frame_index_++;
    //lidar_scans_.push_back(*cloud);
    scan_odoms_.push_back(odom);
    //Mapweighted();
    //SegmentedFrame(cloud);
    
}



void SegmentedMap::Mapaggregate(pcl::PointCloud<PointType>::Ptr cloud0, Eigen::Matrix4f odom){
    
pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>());
     pcl::ApproximateVoxelGrid<PointType> sor1;
     sor1.setInputCloud(cloud0);
     sor1.setLeafSize(0.05f, 0.05f, 0.05f); //0.1f
      sor1.filter(*cloud);
   std::cout<<"cloud size: "<<cloud->size()<<std::endl;
  //pcl::copyPointCloud(*cloud,*current_weighted_cloud_);
  pcl::PointCloud<PointType> seg_cloud = *cloud;
  //SegmentedFrame(seg_cloud.makeShared());
 // lidar_scans_seg_.push_back(seg_cloud);
 // full_lidar_scans_.push_back(*cloud+*ground_cloud);
  //full_cloud_set_.push_back(*cloud);
   for(int i = 0; i < cloud->size(); i++)
   {
       cloud->points[i].intensity = frame_index_;
       
   }
   pcl::PointCloud<PointType>::Ptr ground_cloud(new pcl::PointCloud<PointType>());
   std::vector<int> indices;
    pcl::removeNaNFromPointCloud(*cloud, *ground_cloud, indices);

    std::cout<<"after removenanpoints cloud size: "<<ground_cloud->size()<<std::endl;
    pcl::NormalEstimation<PointType, pcl::Normal> ne;
	ne.setInputCloud (ground_cloud);
	pcl::search::KdTree<PointType>::Ptr tree (new pcl::search::KdTree<PointType> ());
	ne.setSearchMethod (tree);
	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
	ne.setRadiusSearch (0.3);
	//计算特征值
	ne.compute (*cloud_normals);

    for(int i = 0; i < ground_cloud->size();++i)
    {
        Eigen::Vector3d fake_normal = odom.block<3,1>(0,0).cast<double>() - ground_cloud->points.at(i).getVector3fMap().cast<double>();
        
        auto normal_vec = cloud_normals->at(i).getNormalVector3fMap().cast<double>();
        
                if(normal_vec.dot(fake_normal)<0)
                     global_ground_cloud_norm_.push_back(-1 * normal_vec);
                else
                     global_ground_cloud_norm_.push_back(normal_vec);
    }
  //*global_ground_cloud_normals_ += *cloud_normals;
   lidar_ground_scans_.push_back(*ground_cloud);
   *global_ground_map_ += *ground_cloud;

   
    *global_map_ += *cloud;
    frame_index_++;
    //lidar_scans_.push_back(*cloud);
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
      if(segmented_map_.at(i).obs_times_<5)
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
    std::cout<<"global_map_pure_ "<<global_map_pure_->size()<<std::endl;
    //*global_map_pure_ += *global_ground_map_;
    //std::cout<<"ground added global_map_pure_ "<<global_map_pure_->size()<<std::endl;
    std::string n = "/home/wjk/catkin_ws_vbdfusion/results/dynamiccolormap.pcd";
     pcl::io::savePCDFileASCII(n, *rgb_pointcloud); //将点云保存到PCD文件中
     //std::string n1 = "/home/wjk/catkin_ws_vbdfusion/results/globalmap.pcd";
     //pcl::io::savePCDFileASCII(n1, *global_map_pure_); //将点云保存到PCD文件中
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

float xy2theta( const float & _x, const float & _y )
{
    if ( _x >= 0 & _y >= 0) 
        return (180/M_PI) * atan(_y / _x);

    if ( _x < 0 & _y >= 0) 
        return 180 - ( (180/M_PI) * atan(_y / (-_x)) );

    if ( _x < 0 & _y < 0) 
        return 180 + ( (180/M_PI) * atan(_y / _x) );

    if ( _x >= 0 & _y < 0)
        return 360 - ( (180/M_PI) * atan((-_y) / _x) );
} // xy2theta

Eigen::MatrixXd SegmentedMap::makeScancontext( pcl::PointCloud<PointType> & _scan_down )
{
    
     
    int num_pts_scan_down = _scan_down.points.size();

    // main
    const int NO_POINT = -1000;
    Eigen::MatrixXd desc = NO_POINT * Eigen::MatrixXd::Ones(PC_NUM_RING, PC_NUM_SECTOR);
    depth_sc_ = NO_POINT * Eigen::MatrixXd::Ones(PC_NUM_RING, PC_NUM_SECTOR);
    PointType pt;
    float azim_angle, azim_range; // wihtin 2d plane
    int ring_idx, sctor_idx;
    for (int pt_idx = 0; pt_idx < num_pts_scan_down; pt_idx++)
    {
        pt.x = _scan_down.points[pt_idx].x; 
        pt.y = _scan_down.points[pt_idx].y;
        pt.z =  _scan_down.points[pt_idx].z + LIDAR_HEIGHT; // naive adding is ok (all points should be > 0).

        // xyz to ring, sector
        azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);
        azim_angle = xy2theta(pt.x, pt.y);

        // if range is out of roi, pass
        if( azim_range > PC_MAX_RADIUS )
            continue;
        
        ring_idx = std::max( std::min( PC_NUM_RING, int(ceil( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 1 );
        sctor_idx = std::max( std::min( PC_NUM_SECTOR, int(ceil( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 1 );

        // taking maximum z 
        if ( desc(ring_idx-1, sctor_idx-1) < pt.z ) // -1 means cpp starts from 0
            desc(ring_idx-1, sctor_idx-1) = pt.z; // update for taking maximum value at that bin
        
        float depth0  = pt.getVector3fMap().norm();
        if ( depth_sc_(ring_idx-1, sctor_idx-1) < depth0 ) // -1 means cpp starts from 0
            depth_sc_(ring_idx-1, sctor_idx-1) = depth0; // update for taking maximum value at that bin

    }

    // reset no points to zero (for cosine dist later)
    for ( int row_idx = 0; row_idx < desc.rows(); row_idx++ )
        for ( int col_idx = 0; col_idx < desc.cols(); col_idx++ )
        {
            if( desc(row_idx, col_idx) == NO_POINT )
                desc(row_idx, col_idx) = 0;
            if( depth_sc_(row_idx, col_idx) == NO_POINT )
                depth_sc_(row_idx, col_idx) = 0;
        }
            

    return desc;
} // SCManager::makeScancontext


Eigen::MatrixXd SegmentedMap::makeMeanScancontext( pcl::PointCloud<PointType> & _scan_down )
{
    int num_pts_scan_down = _scan_down.points.size();

    // main
    const int NO_POINT = -1000;
    Eigen::MatrixXd desc = 0 * Eigen::MatrixXd::Ones(PC_NUM_RING, PC_NUM_SECTOR);
    Eigen::MatrixXd desc_count = 0 * Eigen::MatrixXd::Ones(PC_NUM_RING, PC_NUM_SECTOR);

    PointType pt;
    float azim_angle, azim_range; // wihtin 2d plane
    int ring_idx, sctor_idx;
    for (int pt_idx = 0; pt_idx < num_pts_scan_down; pt_idx++)
    {
        pt.x = _scan_down.points[pt_idx].x; 
        pt.y = _scan_down.points[pt_idx].y;
        pt.z =  _scan_down.points[pt_idx].z + LIDAR_HEIGHT; // naive adding is ok (all points should be > 0).

        // xyz to ring, sector
        azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);
        azim_angle = xy2theta(pt.x, pt.y);

        // if range is out of roi, pass
        if( azim_range > PC_MAX_RADIUS )
            continue;
        
        ring_idx = std::max( std::min( PC_NUM_RING, int(ceil( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 1 );
        sctor_idx = std::max( std::min( PC_NUM_SECTOR, int(ceil( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 1 );

        // taking maximum z 
        desc(ring_idx-1, sctor_idx-1) += pt.z;
        desc_count(ring_idx-1, sctor_idx-1)++;
        
    }

    // reset no points to zero (for cosine dist later)
    for ( int row_idx = 0; row_idx < desc.rows(); row_idx++ )
        for ( int col_idx = 0; col_idx < desc.cols(); col_idx++ )
            if( desc(row_idx, col_idx) > 0 )
                desc(row_idx, col_idx) = desc(row_idx, col_idx)/desc_count(row_idx, col_idx);
    
    

    return desc;
} // SCManager::makeMeanScancontext


std::vector<std::vector<float>>  SegmentedMap::detectObs(Eigen::MatrixXd scan_sc, Eigen::MatrixXd map_sc)
{
               std::vector<std::vector<float>> dyn_block;
               dyn_block.clear();
               dyn_block.resize(60);
               for(int k = 0; k < scan_sc.cols(); ++k)
              {
                    for(int m = 0; m < scan_sc.rows(); ++m)
                    {
                               float h =   scan_sc(m,k);
                               float h0 = map_sc(m,k);
                               if(h==0)
                                   continue;
                              if(h < 0.2 and h0 > 0.5)
                              {
                                  float depth0 = depth_sc_(m,k);
                                  dyn_block.at(k).push_back(depth0);
                                  break;
                              } 
                    }
                }
                return dyn_block;
}

void SegmentedMap::Getcol(PointType pt,int& row_ind, int& col_ind)
{
   float azim_angle = xy2theta(pt.x, pt.y);
   float azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);
  int ring_idx = std::max( std::min( PC_NUM_RING, int(ceil( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 1 );
  int sctor_idx = std::max( std::min( PC_NUM_SECTOR, int(ceil( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 1 );
  //int id = (ring_idx-1)*PC_NUM_SECTOR+ sctor_idx-1;
  row_ind =  (ring_idx-1);
  col_ind = sctor_idx-1;
}

/*
int SegmentedMap::Getcol(PointType pt)
{
   float azim_angle = xy2theta(pt.x, pt.y);
   float azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);
  int ring_idx = std::max( std::min( PC_NUM_RING, int(ceil( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 1 );
  int sctor_idx = std::max( std::min( PC_NUM_SECTOR, int(ceil( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 1 );
  int id = (ring_idx-1)*PC_NUM_SECTOR+ sctor_idx-1;
  return id;
}*/