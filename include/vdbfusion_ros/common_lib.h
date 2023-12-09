#ifndef COMMON_LIB_H
#define COMMON_LIB_H

#include <Eigen/Eigen>
#include <eigen_conversions/eigen_msg.h>
#include <nav_msgs/Odometry.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <so3_math.h>
#include <tf/transform_broadcaster.h>


using namespace std;
using namespace Eigen;

#define PI_M (3.14159265358)
#define G_m_s2 (9.81)   // Gravaty const in GuangDong/China
#define DIM_STATE (18)  // Dimension of states (Let Dim(SO(3)) = 3)
#define DIM_PROC_N (12) // Dimension of process noise (Let Dim(SO(3)) = 3)
#define CUBE_LEN (6.0)
#define LIDAR_SP_LEN (2)
// old init
#define INIT_COV (0.0000001)
#define NUM_MATCH_POINTS (5)
#define MAX_MEAS_DIM (10000)

#define VEC_FROM_ARRAY(v) v[0], v[1], v[2]
#define MAT_FROM_ARRAY(v) v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]
#define CONSTRAIN(v, min, max) ((v > min) ? ((v < max) ? v : max) : min)
#define ARRAY_FROM_EIGEN(mat) mat.data(), mat.data() + mat.rows() * mat.cols()
#define STD_VEC_FROM_EIGEN(mat)                                                \
  vector<decltype(mat)::Scalar>(mat.data(),                                    \
                                mat.data() + mat.rows() * mat.cols())
#define DEBUG_FILE_DIR(name) (string(string(ROOT_DIR) + "Log/" + name))

typedef pcl::PointXYZINormal PointNType;
typedef pcl::PointCloud<PointType> PointCloudXYZI;
typedef vector<PointType, Eigen::aligned_allocator<PointType>> PointVector;
typedef Vector3d V3D;
typedef Matrix3d M3D;
typedef Vector3f V3F;
typedef Matrix3f M3F;

#define MD(a, b) Matrix<double, (a), (b)>
#define VD(a) Matrix<double, (a), 1>
#define MF(a, b) Matrix<float, (a), (b)>
#define VF(a) Matrix<float, (a), 1>

M3D Eye3d(M3D::Identity());
M3F Eye3f(M3F::Identity());
V3D Zero3d(0, 0, 0);
V3F Zero3f(0, 0, 0);
Vector3d Lidar_offset_to_IMU(0, 0, 0);





template <typename T> T rad2deg(T radians) { return radians * 180.0 / PI_M; }

template <typename T> T deg2rad(T degrees) { return degrees * PI_M / 180.0; }



/* comment
plane equation: Ax + By + Cz + D = 0
convert to: A/D*x + B/D*y + C/D*z = -1
solve: A0*x0 = b0
where A0_i = [x_i, y_i, z_i], x0 = [A/D, B/D, C/D]^T, b0 = [-1, ..., -1]^T
normvec:  normalized x0
*/
template <typename T>
bool esti_normvector(Matrix<T, 3, 1> &normvec, const PointVector &point,
                     const T &threshold, const int &point_num) {
  MatrixXf A(point_num, 3);
  MatrixXf b(point_num, 1);
  b.setOnes();
  b *= -1.0f;

  for (int j = 0; j < point_num; j++) {
    A(j, 0) = point[j].x;
    A(j, 1) = point[j].y;
    A(j, 2) = point[j].z;
  }
  normvec = A.colPivHouseholderQr().solve(b);

  for (int j = 0; j < point_num; j++) {
    if (fabs(normvec(0) * point[j].x + normvec(1) * point[j].y +
             normvec(2) * point[j].z + 1.0f) > threshold) {
      return false;
    }
  }

  normvec.normalize();
  return true;
}

template <typename T>
bool esti_plane(Matrix<T, 4, 1> &pca_result, const PointVector &point,
                const T &threshold) {
  Matrix<T, NUM_MATCH_POINTS, 3> A;
  Matrix<T, NUM_MATCH_POINTS, 1> b;
  b.setOnes();
  b *= -1.0f;

  for (int j = 0; j < NUM_MATCH_POINTS; j++) {
    A(j, 0) = point[j].x;
    A(j, 1) = point[j].y;
    A(j, 2) = point[j].z;
  }

  Matrix<T, 3, 1> normvec = A.colPivHouseholderQr().solve(b);

  for (int j = 0; j < NUM_MATCH_POINTS; j++) {
    if (fabs(normvec(0) * point[j].x + normvec(1) * point[j].y +
             normvec(2) * point[j].z + 1.0f) > threshold) {
      return false;
    }
  }

  T n = normvec.norm();
  pca_result(0) = normvec(0) / n;
  pca_result(1) = normvec(1) / n;
  pca_result(2) = normvec(2) / n;
  pca_result(3) = 1.0 / n;
  return true;
}

#endif