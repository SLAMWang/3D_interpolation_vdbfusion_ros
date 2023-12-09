#ifndef ROADDETECT_H
#define ROADDETECT_H

#include "utility.h"

#include <string>
#include <queue>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_ros/point_cloud.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

class roadDetect{

public:
    roadDetect(const std::string &configPath="default");
    void detectSingleCloud(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud,pcl::PointCloud<pcl::PointXYZI>::Ptr roadCloud);
    void buildRoadGrid(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud,const std::vector<int> &index,gridSubMap &gridSubMap);
    void detectDynamicCloud(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud,std::vector<int> &dynamicIndex,gridSubMap &gridSubMap,std::vector<int> &dynamicGridIndex);
    void setConfig(const std::string &configPath);

private:
    int lineWidth_,scanNumber_,reverseThreshold_,corrosionWidth_,roadPointThreshold_;
    double heightRatio_,lineRatio_,heightThreshold_;



    std::vector<int> grid;
    void flood(std::vector<double> &heights, const int &xSize, const int &ySize);
    void corrosion(const int &xSize, const int &ySize);
};


#endif // ROADDETECT_H
