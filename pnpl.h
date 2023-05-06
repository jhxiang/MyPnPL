#ifndef __PNPL_H__
#define __PNPL_H__
#include <opencv2/opencv.hpp>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_edge.h>
#include <g2o/core/optimizable_graph.h>
#include <g2o/core/batch_stats.h>


// flags: 1:only lines 2:only points 3:points and lines
void MonoPnPL(const std::vector<cv::Vec6f> &lns3d, const std::vector<cv::Vec4f> &lns2d,
            const std::vector<cv::Point3f> &pts3d, const std::vector<cv::Point2f> &pts2d, 
            const cv::Mat &K, cv::Mat &R, std::vector<float> &oula, cv::Mat &t, 
            Eigen::Vector3d trans, Eigen::Quaterniond q, int flags);

void BinoPnPL(const std::vector<cv::Vec6f> &lns3d, const std::vector<cv::Vec4f> &l_lns2d, const std::vector<cv::Vec4f> &r_lns2d,
              const std::vector<cv::Point3f> &pts3d, const std::vector<cv::Point2f> &l_pts2d, const std::vector<cv::Point2f> &r_pts2d,
              const cv::Mat &K_l, const cv::Mat &K_r, const Eigen::Matrix3d R_L2R, const Eigen::Matrix<double, 3, 1> T_L2R,
              cv::Mat &R, std::vector<float> &oula, cv::Mat &t, Eigen::Vector3d trans, Eigen::Quaterniond q, int flags);
#endif
