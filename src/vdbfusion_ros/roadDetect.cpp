#include "roadDetect.h"

roadDetect::roadDetect(const std::string &configPath){
    std::string lidarPath, sequence, calibPath, gtPath;

    if (configPath=="default"){
        lineWidth_=10;
        scanNumber_=3;
        lineRatio_=0.33;
        heightRatio_=0.66;
        reverseThreshold_=30;
        corrosionWidth_=2;
        roadPointThreshold_=10;
        heightThreshold_=0.2;
    }
    else{
        boost::property_tree::ptree pt;
        boost::property_tree::read_xml(configPath, pt);
        pt = pt.get_child("config");

        lineWidth_=pt.get<int>("lineWidth",10);
        scanNumber_=pt.get<int>("scanNumber",3);
        lineRatio_=pt.get<double>("lineRatio",0.33);
        heightRatio_=pt.get<double>("heightRatio",0.66);
        reverseThreshold_=pt.get<int>("reverseThreshold",30);
        corrosionWidth_=pt.get<int>("corrosionWidth",2);
        roadPointThreshold_=pt.get<int>("roadPointThreshold",10);
        heightThreshold_=pt.get<double>("heightThreshold",0.2);
    }


}

void roadDetect::setConfig(const std::string &configPath){
    boost::property_tree::ptree pt;
    boost::property_tree::read_xml(configPath, pt);
    pt = pt.get_child("config");

    lineWidth_=pt.get<int>("lineWidth",10);
    scanNumber_=pt.get<int>("scanNumber",3);
    lineRatio_=pt.get<double>("lineRatio",0.33);
    heightRatio_=pt.get<double>("heightRatio",0.66);
    reverseThreshold_=pt.get<int>("reverseThreshold",30);
    corrosionWidth_=pt.get<int>("corrosionWidth",2);
    roadPointThreshold_=pt.get<int>("roadPointThreshold",10);
    heightThreshold_=pt.get<double>("heightThreshold",0.2);
}

void roadDetect::detectSingleCloud(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, pcl::PointCloud<pcl::PointXYZI>::Ptr roadCloud){
    roadCloud->resize(0);
    int scanSize(1920);

    std::vector<int> index(16 * scanSize, -1), order{-15, 1, -13, 3, -11, 5, -9, 7, -7, 9, -5, 11, -3, 13, -1, 15},
    pointIndex(cloud->size(), -1);
    int orderNum = 0, indexNow = 0;

    for (int i = 0; i < cloud->size(); i++) {
        auto point = cloud->at(i);
        double theta, radius;
        radius = sqrt(point.x * point.x + point.y * point.y);
        theta = atan2(point.z, radius) * 180 / M_PI;
        while (std::abs(theta - order[orderNum]) > 0.01) {
            indexNow++;
            orderNum = (orderNum + 1) % 16;
        }
        index[indexNow] = i;
        pointIndex[i] = indexNow;
        indexNow++;
        orderNum = (orderNum + 1) % 16;
    }

    std::vector<double> distanceData(scanNumber_ * scanSize, 0), lineDistanceData, heightData;
    for (int scanId = 0; scanId < 2 * scanNumber_; scanId += 2)
        for (int i = 0; i < scanSize; i++) {
            int handleNumber = i + scanId * scanSize / 2;
            if (index[i * 16 + scanId] >= 0) {
                auto point = cloud->at(index[i * 16 + scanId]);
                auto poseI = point.getVector3fMap();
                heightData.push_back(point.z);
                for (auto j = -lineWidth_; j < lineWidth_; j++) {
                    int neighbourNum = 16 * ((i + j + scanSize) % scanSize) + scanId;
                    if (index[neighbourNum] >= 0) {
                        auto poseN = cloud->at(index[neighbourNum]).getVector3fMap();
                        distanceData[handleNumber] += (poseI - poseN).norm();
                    }
                }
                if (point.z < -0.5)
                    lineDistanceData.push_back(distanceData[handleNumber]);
            }
        }


    auto m = lineDistanceData.begin() + static_cast<int>(lineDistanceData.size() * lineRatio_);
    std::nth_element(lineDistanceData.begin(), m, lineDistanceData.end());
    m = heightData.begin() + static_cast<int>(heightData.size() * heightRatio_);
    std::nth_element(heightData.begin(), m, heightData.end());
    double heightThreshold, lineDistanceThreshold;
    heightThreshold = heightData[static_cast<int>(heightData.size() * heightRatio_)];
    lineDistanceThreshold = lineDistanceData[static_cast<int>(lineDistanceData.size() * lineRatio_)];

    for (int scanId = 0; scanId < 2 * scanNumber_; scanId += 2)
        for (int i = 0; i < scanSize; i++)
            if (index[i * 16 + scanId] >= 0) {
                int handleNumber = i + scanId * scanSize / 2;
                auto point = cloud->at(index[i * 16 + scanId]);
                if (distanceData[handleNumber] < lineDistanceThreshold && point.z < heightThreshold) {
                    point.intensity=distanceData[handleNumber];
                    roadCloud->push_back(point);
                }
            }

}

void roadDetect::flood(std::vector<double> &heights, const int &xSize, const int &ySize) {

    for (auto &iter: grid)
        if (iter > roadPointThreshold_)
            iter = 1;
        else
            iter = 0;

    std::vector<std::vector<int>> neighbor{{1,  0},{-1, 0},{0,  -1},{0,  1}};

    std::vector<int> visit(grid.size(), 0);
    for (int i = 0; i < grid.size(); i++) {
        if (!visit[i]) {
            visit[i]=1;
            std::vector<int> clusterList;
            std::queue<std::pair<int,int>> searchList;
            searchList.push(std::pair<int,int>(i / ySize,i % ySize));
            while (!searchList.empty()) {
                std::pair<int,int> handlePoint=searchList.front();
                searchList.pop();
                int xx,yy,pos;
                xx=handlePoint.first;
                yy=handlePoint.second;
                pos=xx*ySize+yy;
                clusterList.push_back(pos);
                for (int i = 0; i < 4; i++) {
                    int x, y;
                    x = xx + neighbor[i][0];
                    y = yy + neighbor[i][1];
                    if (x >= 0 && x < xSize && y >= 0 && y < ySize && visit[x * ySize + y]==0 && grid[pos] == grid[x * ySize + y]){
                        visit[x * ySize + y]=1;
                        searchList.push(std::pair<int,int>(x,y));
                    }
                }
            }
            if (clusterList.size() < reverseThreshold_) {
                double averageHeight = 0;
                int nums = 0;
                for (auto iter: clusterList) {
                    if (!grid[iter]) {
                        int xx, yy;
                        xx = iter / ySize;
                        yy = iter % ySize;
                        for (int x = -1; x < 2; x++)
                            for (int y = -1; y < 2; y++)
                                if (xx + x >= 0 && xx + x < xSize && yy + y >= 0 && yy + y < ySize &&
                                        grid[(xx + x) * ySize + yy + y]) {
                                    averageHeight += heights[(xx + x) * ySize + yy + y];
                                    nums++;
                                }
                    }
                    grid[iter] = 1 - grid[iter];
                }
                averageHeight /= nums;
                if (grid[i])
                    for (auto iter: clusterList) {
                        heights[iter] = averageHeight;
                    }
            }
        }

    }
}

void roadDetect::corrosion(const int &xSize, const int &ySize) {

    for (int i = 0; i < grid.size(); i++)
        if (grid[i]) {
            int xx, yy;
            xx = i / ySize;
            yy = i % ySize;
            int zeroNum(0);
            for (int x = -corrosionWidth_; x < corrosionWidth_+1; x++)
                for (int y = -corrosionWidth_; y < corrosionWidth_+1; y++)
                    if (xx + x >= 0 && xx + x < xSize && yy + y >= 0 && yy + y < ySize &&
                            !grid[(xx + x) * ySize + yy + y])
                        zeroNum++;
            if (zeroNum > corrosionWidth_*corrosionWidth_)
                grid[i] = 2;
        }

    for (int i = 0; i < grid.size(); i++)
        if (grid[i] == 2)
            grid[i] = 0;
}

void roadDetect::buildRoadGrid(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud,const std::vector<int> &index,gridSubMap &gridSubMap){

    pcl::PointCloud<pcl::PointXYZI>::Ptr subCloud(new pcl::PointCloud<pcl::PointXYZI>);
    for (int i = 0; i < index.size(); ++i)
        subCloud->push_back(cloud->at(index[i]));

    pcl::PointXYZI min3D,max3D;
    pcl::getMinMax3D(*subCloud,min3D,max3D);
    double xZero = min3D.x, yZero=min3D.y;
    int xSize = (std::ceil(max3D.x/5)*5-xZero) / gridSubMap.resolution, ySize = (std::ceil(max3D.y/5)*5-yZero) / gridSubMap.resolution;

    grid.resize(xSize * ySize, 0);
    std::vector<double> gridAverageHeight(xSize * ySize, 0);
    Eigen::Vector3f centerPoint;
    centerPoint=0.5*(min3D.getVector3fMap()+max3D.getVector3fMap());
    gridSubMap.centerPoint=centerPoint;

    for (auto iter: subCloud->points) {
        int xx, yy;
        xx = std::ceil((iter.x - xZero) / gridSubMap.resolution);
        yy = std::ceil((iter.y - yZero) / gridSubMap.resolution);
        if (xx >= 0 && xx < xSize && yy >= 0 && yy < ySize) {
            grid[xx * ySize + yy]++;
            gridAverageHeight[xx * ySize + yy] += iter.z;
        }
    }


    for (int i = 0; i < xSize * ySize; i++) {
        gridAverageHeight[i] /= grid[i];
    }

    flood(gridAverageHeight, xSize, ySize);
    corrosion(xSize, ySize);
    gridSubMap.gridMap->sensor_origin_ = Eigen::Vector4f(min3D.x, min3D.y,min3D.z, 1);
    gridSubMap.gridMap->sensor_orientation_.x()=gridSubMap.resolution;
    gridSubMap.gridMap->resize(xSize*ySize);
    gridSubMap.gridMap->height = ySize;
    gridSubMap.gridMap->width = xSize;

    gridSubMap.height=ySize;
    gridSubMap.width=xSize;
    gridSubMap.origin=min3D.getArray3fMap();

    for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++)
            if (grid[i*ySize+j]>0)
                gridSubMap.gridMap->at(i, j).intensity = gridAverageHeight[i * ySize + j];
            else
                gridSubMap.gridMap->at(i, j).intensity = -100;
    }

}

void roadDetect::detectDynamicCloud(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud,std::vector<int> &dynamicIndex,gridSubMap &gridSubMap,std::vector<int> &dynamicGridIndex){

    dynamicGridIndex.resize(gridSubMap.width*gridSubMap.height,0);
    double xMin=gridSubMap.origin(0),yMin=gridSubMap.origin(1);
    int subMapSize(0),height(0),occupy(0),disError(0);
    for (int i=0;i<cloud->size();i++) {
        auto iter=cloud->at(i);
        Eigen::Vector3f pose = iter.getVector3fMap();
        auto dis=(pose-gridSubMap.centerPoint).head(2).norm();

        int xx, yy;
        xx = std::ceil((iter.x - xMin) / gridSubMap.resolution);
        yy = std::ceil((iter.y - yMin) / gridSubMap.resolution);

        double heightThreshold = std::max(heightThreshold_, heightThreshold_ / 10 * dis);

        if (xx >= 0 && xx < gridSubMap.width && yy >= 0 && yy < gridSubMap.height && iter.z < gridSubMap.gridMap->at(xx,yy).intensity+2
                &&gridSubMap.gridMap->at(xx,yy).intensity>-50 && iter.z > gridSubMap.gridMap->at(xx,yy).intensity + heightThreshold_) {
            dynamicIndex.push_back(i);
            dynamicGridIndex[xx*gridSubMap.height+yy]=1;
        }
    }

}
