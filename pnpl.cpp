#include <Eigen/Dense>
#include <vector>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/solver.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/types/sba/types_six_dof_expmap.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/types/sba/g2o_types_sba_api.h>

#include "pnpl.h"

#define PI 3.1415926536

using namespace Eigen;
using namespace std;

g2o::Vector2 project(const g2o::Vector3 &v)
{
    g2o::Vector2 res;
    res(0) = v(0) / v(2);
    res(1) = v(1) / v(2);
    return res;
}

g2o::Vector2 mycam_map(const g2o::Vector3 &trans_xyz, const g2o::Vector2 &principle_point, const g2o::Vector2 &focal_length)
{
    g2o::Vector2 proj = project(trans_xyz);
    g2o::Vector2 res;
    res[0] = proj[0] * focal_length[0] + principle_point[0];
    res[1] = proj[1] * focal_length[1] + principle_point[1];
    return res;
}

namespace g2o {
    typedef Eigen::Matrix<double, 8, 1, Eigen::ColMajor> Vector8d;
    typedef Eigen::Matrix<double, 4, 1, Eigen::ColMajor> Vector4d;
    typedef Eigen::Matrix<double, 6, 1, Eigen::ColMajor>  Vector6D;
    class VertexSBALine : public BaseVertex<6, Vector6D>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        VertexSBALine();
        virtual bool read(std::istream& is);
        virtual bool write(std::ostream& os) const;

        virtual void setToOriginImpl() {
            _estimate.fill(0.);
        }

        virtual void oplusImpl(const double* update)
        {
            Eigen::Map<const Vector6D> v(update);
            _estimate += v;
        }
    };

    VertexSBALine::VertexSBALine() : BaseVertex<6, Vector6D>()
    {
    }

    bool VertexSBALine::read(std::istream& is)
    {
        Vector6D lv;
        for (int i = 0; i < 6; i++)
            is >> _estimate[i];
        return true;
    }

    bool VertexSBALine::write(std::ostream& os) const
    {
        Vector6D lv = estimate();
        for (int i = 0; i < 6; i++) {
            os << lv[i] << " ";
        }
        return os.good();
    }

    // ---------------------------------------------一元边，单目优化位姿------------------------------------------
    class EdgeLine : public BaseUnaryEdge<4, Vector4d, VertexSE3Expmap>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        EdgeLine(){}

        bool read(std::istream &is);

        bool write(std::ostream &os) const;

        void computeError()
        {
            const VertexSE3Expmap *v1 = static_cast<const VertexSE3Expmap *>(_vertices[0]);
            // const CameraParameters *cam = static_cast<const CameraParameters *>(parameter(0));
            Vector4d obs(_measurement);
            Vector3d Xc1(v1->estimate().map(Vector3d(_line[0], _line[1], _line[2])));
            Vector3d Xc2(v1->estimate().map(Vector3d(_line[3], _line[4], _line[5])));
            Vector2d u1 = mycam_map(Xc1, _principle_point, _focal_length);
            Vector2d u2 = mycam_map(Xc2, _principle_point, _focal_length);
            double dx = obs[2] - obs[0];
            double dy = obs[3] - obs[1];
            double n = hypot(dx, dy);
            dx /= n;
            dy /= n;
            double d = -dy * obs[0] + dx * obs[1];
            double dist1 = -dy * u1[0] + dx * u1[1] - d;
            double dist2 = -dy * u2[0] + dx * u2[1] - d;
            _error = Vector4d(dist1, dist2, 0, 0);
        }

        void setparams(const Vector6D &line, const Vector2d &focal_length, const Vector2d &principle_point) { 
            _line = line;
            _focal_length = focal_length;
            _principle_point = principle_point;
        }
        // CameraParameters *_cam;

        Vector6D _line;
        Vector2d _focal_length;
        Vector2d _principle_point;
    };

    // EdgeLine::EdgeLine() : BaseUnaryEdge<4, Vector4d, VertexSE3Expmap>()
    // {
    //     _cam = 0;
    //     resizeParameters(1);
    //     installParameter(_cam, 0);
    // }

    bool EdgeLine::read(std::istream &is)
    {
        // int paramId;
        // is >> paramId;
        // setParameterId(0, paramId);

        for (int i = 0; i < 4; i++)
        {
            is >> _measurement[i];
        }
        for (int i = 0; i < 4; i++)
            for (int j = i; j < 4; j++)
            {
                is >> information()(i, j);
                if (i != j)
                    information()(j, i) = information()(i, j);
            }
        return true;
    }

    bool EdgeLine::write(std::ostream &os) const
    {
        // os << _cam->id() << " ";
        for (int i = 0; i < 4; i++)
        {
            os << measurement()[i] << " ";
        }

        for (int i = 0; i < 4; i++)
            for (int j = i; j < 4; j++)
            {
                os << " " << information()(i, j);
            }
        return os.good();
        return true;
    }

    // -------------------------------------------一元边，双目优化位姿(线)-----------------------------------------
    class BinoEdgeLine : public BaseUnaryEdge<8, Vector8d, VertexSE3Expmap>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        BinoEdgeLine(){}

        bool read(std::istream &is);

        bool write(std::ostream &os) const;

        void computeError()
        {
            const VertexSE3Expmap *v1 = static_cast<const VertexSE3Expmap *>(_vertices[0]);
            Vector8d obs(_measurement);
            Vector3d Xc1_l(v1->estimate().map(Vector3d(_line[0], _line[1], _line[2])));
            Vector3d Xc2_l(v1->estimate().map(Vector3d(_line[3], _line[4], _line[5])));
            Vector2d u1 = mycam_map(Xc1_l, _principle_point_l, _focal_length_l);
            Vector2d u2 = mycam_map(Xc2_l, _principle_point_l, _focal_length_l);
            double dx1 = obs[2] - obs[0];
            double dy1 = obs[3] - obs[1];
            double n1 = hypot(dx1, dy1);
            dx1 /= n1;
            dy1 /= n1;
            double d1 = -dy1 * obs[0] + dx1 * obs[1];
            double dist1 = -dy1 * u1[0] + dx1 * u1[1] - d1;
            double dist2 = -dy1 * u2[0] + dx1 * u2[1] - d1;
            
            // 左相机转换为右相机
            Vector3d Xc1_r = R_L2R * Xc1_l + T_L2R;
            Vector3d Xc2_r = R_L2R * Xc2_l + T_L2R;
            Vector2d u3 = mycam_map(Xc1_r, _principle_point_r, _focal_length_r);
            Vector2d u4 = mycam_map(Xc2_r, _principle_point_r, _focal_length_r);
            double dx2 = obs[6] - obs[4];
            double dy2 = obs[7] - obs[5];
            double n2 = hypot(dx2, dy2);
            dx2 /= n2;
            dy2 /= n2;
            double d2 = -dy2 * obs[4] + dx2 * obs[5];
            double dist3 = -dy2 * u3[0] + dx2 * u3[1] - d2;
            double dist4 = -dy2 * u4[0] + dx2 * u4[1] - d2;
            _error = Vector8d(dist1, dist2, 0, 0, dist3, dist4, 0, 0);
        }

        void setparams(const Vector6D &line, const Vector2d &focal_length_l, const Vector2d &focal_length_r,
                       const Vector2d &principle_point_l, const Vector2d &principle_point_r) { 
            _line = line;
            _focal_length_l = focal_length_l;
            _focal_length_r = focal_length_r;
            _principle_point_l = principle_point_l;
            _principle_point_r = principle_point_r;
        }
        // CameraParameters *_cam;

        Vector6D _line;
        Matrix<double, 3, 3> R_L2R; // 左相机到右相机的变换矩阵
        Matrix<double, 3, 1> T_L2R;
        Vector2d _focal_length_l;
        Vector2d _principle_point_l;
        Vector2d _focal_length_r;
        Vector2d _principle_point_r;
    };

    // BinoEdgeLine::BinoEdgeLine() : BaseUnaryEdge<8, Vector8d, VertexSE3Expmap>()
    // {
    //     _cam = 0;
    //     resizeParameters(1);
    //     installParameter(_cam, 0);
    // }

    bool BinoEdgeLine::read(std::istream &is)
    {
        // int paramId;
        // is >> paramId;
        // setParameterId(0, paramId);

        for (int i = 0; i < 8; i++)
        {
            is >> _measurement[i];
        }
        for (int i = 0; i < 8; i++)
            for (int j = i; j < 8; j++)
            {
                is >> information()(i, j);
                if (i != j)
                    information()(j, i) = information()(i, j);
            }
        return true;
    }

    bool BinoEdgeLine::write(std::ostream &os) const
    {
        // os << _cam->id() << " ";
        for (int i = 0; i < 8; i++)
        {
            os << measurement()[i] << " ";
        }

        for (int i = 0; i < 8; i++)
            for (int j = i; j < 8; j++)
            {
                os << " " << information()(i, j);
            }
        return os.good();
        return true;
    }
    
    // ---------------------------------------一元边，双目优化位姿(点)-----------------------------------------
    class EdgeSteoSE3ProjectXYZOnlyPose: public BaseUnaryEdge<4, Vector4d, VertexSE3Expmap>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        EdgeSteoSE3ProjectXYZOnlyPose(){}
        bool read(std::istream &is);
        bool write(std::ostream &os) const;
        void computeError();

    public:
        Vector2d cam_project_l(const Vector3d &trans_xyz) const;
        Vector2d cam_project_r(const Vector3d &trans_xyz) const;

        Vector3d Xw;
        Matrix<double, 3, 3> R_L2R; // 左相机到右相机的变换矩阵
        Matrix<double, 3, 1> T_L2R;
        double fx_l, fy_l, cx_l, cy_l, fx_r, fy_r, cx_r, cy_r;
    };

    bool EdgeSteoSE3ProjectXYZOnlyPose::read(std::istream &is)
    {
        internal::readVector(is, _measurement);
        return readInformationMatrix(is);
    }

    bool EdgeSteoSE3ProjectXYZOnlyPose::write(std::ostream &os) const
    {
        internal::writeVector(os, measurement());
        return writeInformationMatrix(os);
    }

    Vector2d EdgeSteoSE3ProjectXYZOnlyPose::cam_project_l(const Vector3d &trans_xyz) const
    {
        Vector2d proj = project(trans_xyz);
        Vector2d res;
        res[0] = proj[0] * fx_l + cx_l;
        res[1] = proj[1] * fy_l + cy_l;
        return res;
    }

    Vector2d EdgeSteoSE3ProjectXYZOnlyPose::cam_project_r(const Vector3d &trans_xyz) const
    {
        Vector3d trans_xyz_r = R_L2R * trans_xyz + T_L2R;
        Vector2d proj = project(trans_xyz_r);
        Vector2d res;
        res[0] = proj[0] * fx_r + cx_r;
        res[1] = proj[1] * fy_r + cy_r;
        return res;
    }

    void EdgeSteoSE3ProjectXYZOnlyPose::computeError()
    {
        const VertexSE3Expmap *v1 = static_cast<const VertexSE3Expmap *>(_vertices[0]);
        Vector4d obs(_measurement);
        Vector2d error1 = cam_project_l(v1->estimate().map(Xw));
        Vector2d error2 = cam_project_r(v1->estimate().map(Xw));
        _error = obs - Vector4d(error1[0], error1[1], error2[0], error2[1]);
    }


    //----------------------------------------二元边，优化位姿和三维模型---------------------------------------------
    class EdgeProjectLine : public  BaseBinaryEdge<4, Vector4d, VertexSBALine, VertexSE3Expmap> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        EdgeProjectLine();

        bool read(std::istream& is);

        bool write(std::ostream& os) const;

        void computeError() {
            const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
            const VertexSBALine* v2 = static_cast<const VertexSBALine*>(_vertices[0]);
            const CameraParameters* cam
                = static_cast<const CameraParameters*>(parameter(0));
            Vector4d obs(_measurement);
            Vector6D est = v2->estimate();
            Vector2d u1 = cam->cam_map(v1->estimate().map(Vector3d(est[0], est[1], est[2])));
            Vector2d u2 = cam->cam_map(v1->estimate().map(Vector3d(est[3], est[4], est[5])));
            double dx = obs[2] - obs[0];
            double dy = obs[3] - obs[1];
            double n = hypot(dx, dy);
            dx /= n;
            dy /= n;
            double d = -dy * obs[0] + dx * obs[1];
            double dist1 = -dy * u1[0] + dx * u1[1] - d;
            double dist2 = -dy * u2[0] + dx * u2[1] - d;
            _error = Vector4d(dist1, dist2, 0, 0);
        }

        CameraParameters* _cam;
    };

    EdgeProjectLine::EdgeProjectLine() : BaseBinaryEdge<4, Vector4d, VertexSBALine, VertexSE3Expmap>() {
        _cam = 0;
        resizeParameters(1);
        installParameter(_cam, 0);
    }

    bool EdgeProjectLine::read(std::istream& is) {
        int paramId;
        is >> paramId;
        setParameterId(0, paramId);

        for (int i = 0; i < 4; i++) {
            is >> _measurement[i];
        }
        for (int i = 0; i < 4; i++)
            for (int j = i; j < 4; j++) {
                is >> information()(i, j);
                if (i != j)
                    information()(j, i) = information()(i, j);
            }
        return true;
    }

    bool EdgeProjectLine::write(std::ostream& os) const {
        os << _cam->id() << " ";
        for (int i = 0; i < 4; i++) {
            os << measurement()[i] << " ";
        }

        for (int i = 0; i < 4; i++)
            for (int j = i; j < 4; j++) {
                os << " " << information()(i, j);
            }
        return os.good();
        return true;
    }
} //namespace g2o

vector<double> R2angle(cv::Mat_<double> R)
    {
        vector<double> angle;
        double angle1, angle2, angle3;

        // 1:俯仰 2:方位 3:滚转
        angle3 = atan(R(0, 1) / R(1, 1));
        angle1 = asin(-R(2, 1));
        angle2 = asin(-R(2, 0) / cos(angle1));
        angle.push_back(angle1 * 180 / PI);
        angle.push_back(angle2 * 180 / PI);
        angle.push_back(angle3 * 180 / PI);
        return angle;
    }

   
void MonoPnPL(const std::vector<cv::Vec6f> &lns3d, const std::vector<cv::Vec4f> &lns2d,
            const std::vector<cv::Point3f> &pts3d, const std::vector<cv::Point2f> &pts2d,
            const cv::Mat &K, cv::Mat &R, std::vector<double> &oula, cv::Mat &t, 
            Eigen::Vector3d trans, Eigen::Quaterniond q, int flags)
{
    int nlns = lns3d.size();
    int npts = pts3d.size();
    // init g2o
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<6, 1>> Block;
    Block::LinearSolverType *linearSolver;
    linearSolver = new g2o::LinearSolverDense<Block::PoseMatrixType>();
    Block *solver_ptr = new Block(unique_ptr<Block::LinearSolverType>(linearSolver));
    g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(unique_ptr<Block>(solver_ptr));
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(false);
    optimizer.setAlgorithm(solver);

    // add vertex for pose
    int vertex_id = 0;
    g2o::SE3Quat pose(q, trans);
    g2o::VertexSE3Expmap *v_se3 = new g2o::VertexSE3Expmap();
    v_se3->setId(vertex_id++);
    v_se3->setEstimate(pose);
    optimizer.addVertex(v_se3);

    // set camera intrinsic
    cv::Mat _K;
    K.convertTo(_K, CV_64FC1);
    // g2o::CameraParameters *cam_params = new g2o::CameraParameters(_K.at<double>(0, 0), Eigen::Vector2d(_K.at<double>(0, 2), _K.at<double>(1, 2)), 0.);
    // cam_params->setId(0);
    // optimizer.addParameter(cam_params);
    
    switch(flags){
        case 1:
            cout << "mono: optimization by lines only" << endl;
            // add vertex for lines and edges for projection
            for (int i = 0; i < nlns; i++)
            {
                const auto &ln3d = lns3d[i];
                const auto &ln2d = lns2d[i];

                g2o::Vector6D temp;
                Eigen::Vector2d focal_length(_K.at<double>(0, 0), _K.at<double>(1, 1));
                Eigen::Vector2d principle_point(_K.at<double>(0, 2), _K.at<double>(1, 2));
                temp << ln3d[0], ln3d[1], ln3d[2], ln3d[3], ln3d[4], ln3d[5];

                g2o::EdgeLine *e = new g2o::EdgeLine();
                e->setparams(temp, focal_length, principle_point);
                e->setVertex(0, v_se3);
                e->setMeasurement(Eigen::Vector4d(ln2d[0], ln2d[1], ln2d[2], ln2d[3]));
                e->information() = Eigen::Matrix4d::Identity();
                e->setRobustKernel(new g2o::RobustKernelHuber);
                // e->setParameterId(0, 0);
                optimizer.addEdge(e);
            }
            break;

        case 2:
            cout << "mono: optimization by points only" << endl;
            // add add vertex for points and edges for projection
            for (int i = 0; i < npts; i++)
            {
                const auto &pt = pts3d[i];
                const auto &uv = pts2d[i];

                g2o::EdgeSE3ProjectXYZOnlyPose *e = new g2o::EdgeSE3ProjectXYZOnlyPose();
                e->setVertex(0, v_se3);
                e->setMeasurement(Eigen::Vector2d(uv.x, uv.y));
                e->information() = Eigen::Matrix2d::Identity();
                e->setRobustKernel(new g2o::RobustKernelHuber);
                // e->setParameterId(0, 0);

                e->fx = _K.at<double>(0, 0);
                e->fy = _K.at<double>(1, 1);
                e->cx = _K.at<double>(0, 2);
                e->cy = _K.at<double>(1, 2);
                e->Xw[0] = pt.x;
                e->Xw[1] = pt.y;
                e->Xw[2] = pt.z;
                optimizer.addEdge(e);
            }
            break;

        case 3:
            cout << "mono: optimization by lines and points" << endl;
            // add add vertex for points and edges for projection
            for (int i = 0; i < npts; i++)
            {
                const auto &pt = pts3d[i];
                const auto &uv = pts2d[i];

                g2o::EdgeSE3ProjectXYZOnlyPose *e = new g2o::EdgeSE3ProjectXYZOnlyPose();
                e->setVertex(0, v_se3);
                e->setMeasurement(Eigen::Vector2d(uv.x, uv.y));
                e->information() = Eigen::Matrix2d::Identity();
                e->setRobustKernel(new g2o::RobustKernelHuber);
                // e->setParameterId(0, 0);

                e->fx = _K.at<double>(0, 0);
                e->fy = _K.at<double>(1, 1);
                e->cx = _K.at<double>(0, 2);
                e->cy = _K.at<double>(1, 2);
                e->Xw[0] = pt.x;
                e->Xw[1] = pt.y;
                e->Xw[2] = pt.z;
                optimizer.addEdge(e);
            }

            // add vertex for lines and edges for projection
            for (int i = 0; i < nlns; i++)
            {
                const auto &ln3d = lns3d[i];
                const auto &ln2d = lns2d[i];

                g2o::Vector6D temp;
                Eigen::Vector2d focal_length(_K.at<double>(0, 0), _K.at<double>(1, 1));
                Eigen::Vector2d principle_point(_K.at<double>(0, 2), _K.at<double>(1, 2));
                temp << ln3d[0], ln3d[1], ln3d[2], ln3d[3], ln3d[4], ln3d[5];

                g2o::EdgeLine *e = new g2o::EdgeLine();
                e->setparams(temp, focal_length, principle_point);
                e->setVertex(0, v_se3);
                e->setMeasurement(Eigen::Vector4d(ln2d[0], ln2d[1], ln2d[2], ln2d[3]));
                e->information() = Eigen::Matrix4d::Identity();
                e->setRobustKernel(new g2o::RobustKernelHuber);
                // e->setParameterId(0, 0);
                optimizer.addEdge(e);
            }
            break;
        
        default:
            cerr << "wrong flags!!!  please input 1 or 2 or 3!" << endl;
            break;
    }
    // optimize
    optimizer.initializeOptimization();
    optimizer.optimize(20);
    // output result
    Eigen::MatrixXd T = Eigen::Isometry3d(v_se3->estimate()).matrix();
    R = (cv::Mat_<double>(3, 3) << T(0, 0), T(0, 1), T(0, 2),
                                   T(1, 0), T(1, 1), T(1, 2),
                                   T(2, 0), T(2, 1), T(2, 2));           
    oula = R2angle(R);
    t = (cv::Mat_<double>(3, 1) << T(0, 3), T(1, 3), T(2, 3));
    t = R.t() * t;
}

void BinoPnPL(const std::vector<cv::Vec6f> &lns3d, const std::vector<cv::Vec4f> &l_lns2d, const std::vector<cv::Vec4f> &r_lns2d,
              const std::vector<cv::Point3f> &pts3d, const std::vector<cv::Point2f> &l_pts2d, const std::vector<cv::Point2f> &r_pts2d,
              const cv::Mat &K_l, const cv::Mat &K_r, const Eigen::Matrix3d R_L2R, const Eigen::Matrix<double, 3, 1> T_L2R, cv::Mat &R, 
              std::vector<double> &oula, cv::Mat &t, Eigen::Vector3d trans, Eigen::Quaterniond q, int flags)
{
    int nlns = lns3d.size();
    int npts = pts3d.size();
    // init g2o
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<6, 1>> Block;
    Block::LinearSolverType *linearSolver;
    linearSolver = new g2o::LinearSolverDense<Block::PoseMatrixType>();
    Block *solver_ptr = new Block(unique_ptr<Block::LinearSolverType>(linearSolver));
    g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(unique_ptr<Block>(solver_ptr));
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(false);
    optimizer.setAlgorithm(solver);

    // add vertex for pose
    int vertex_id = 0;
    g2o::SE3Quat pose(q, trans);
    g2o::VertexSE3Expmap *v_se3 = new g2o::VertexSE3Expmap();
    v_se3->setId(vertex_id++);
    v_se3->setEstimate(pose);
    optimizer.addVertex(v_se3);

    // set camera intrinsic，转换到哪个相机下就用哪个相机的内参
    cv::Mat _K_l, _K_r;
    K_l.convertTo(_K_l, CV_64FC1);
    K_r.convertTo(_K_r, CV_64FC1);
    // g2o::CameraParameters *cam_params = new g2o::CameraParameters(_K_l.at<double>(0, 0), Eigen::Vector2d(_K_l.at<double>(0, 2), _K_l.at<double>(1, 2)), 0.);
    // cam_params->setId(0);
    // optimizer.addParameter(cam_params);

    switch (flags)
    {
    case 1:
            cout << "bino: optimization by lines only" << endl;
            // add vertex for lines and edges for projection
            for (int i = 0; i < nlns; i++)
            {
                const auto &ln3d = lns3d[i];
                const auto &l_ln2d = l_lns2d[i];
                const auto &r_ln2d = r_lns2d[i];

                g2o::Vector6D temp;
                Vector2d focal_length_l(_K_l.at<double>(0, 0), _K_l.at<double>(1, 1));
                Vector2d principle_point_l(_K_l.at<double>(0, 2), _K_l.at<double>(1, 2));
                Vector2d focal_length_r(_K_r.at<double>(0, 0), _K_r.at<double>(1, 1));
                Vector2d principle_point_r(_K_r.at<double>(0, 2), _K_r.at<double>(1, 2));
                temp << ln3d[0], ln3d[1], ln3d[2], ln3d[3], ln3d[4], ln3d[5];

                g2o::BinoEdgeLine *e = new g2o::BinoEdgeLine();
                e->setparams(temp, focal_length_l, focal_length_r, principle_point_l, principle_point_r);
                e->setVertex(0, v_se3);
                e->R_L2R = R_L2R;
                e->T_L2R = T_L2R;
                e->setMeasurement(g2o::Vector8d(l_ln2d[0], l_ln2d[1], l_ln2d[2], l_ln2d[3], r_ln2d[0], r_ln2d[1], r_ln2d[2], r_ln2d[3]));
                e->information() = Eigen::Matrix<double, 8, 8>::Identity();
                e->setRobustKernel(new g2o::RobustKernelHuber);
                // e->setParameterId(0, 0);
                optimizer.addEdge(e);
            }
            break;

    case 2:
            cout << "bino: optimization by points only" << endl;
            // add add vertex for points and edges for projection
            for (int i = 0; i < npts; i++)
            {
                const auto &pt = pts3d[i];
                const auto &l_uv = l_pts2d[i];
                const auto &r_uv = r_pts2d[i];

                g2o::EdgeSteoSE3ProjectXYZOnlyPose *e = new g2o::EdgeSteoSE3ProjectXYZOnlyPose();
                e->setVertex(0, v_se3);
                e->setMeasurement(Eigen::Vector4d(l_uv.x, l_uv.y, r_uv.x, r_uv.y));
                e->information() = Eigen::Matrix4d::Identity();
                // e->setParameterId(0, 0);
                e->setRobustKernel(new g2o::RobustKernelHuber);

                e->fx_l = _K_l.at<double>(0, 0);
                e->fy_l = _K_l.at<double>(1, 1);
                e->cx_l = _K_l.at<double>(0, 2);
                e->cy_l = _K_l.at<double>(1, 2);
                e->fx_r = _K_r.at<double>(0, 0);
                e->fy_r = _K_r.at<double>(1, 1);
                e->cx_r = _K_r.at<double>(0, 2);
                e->cy_r = _K_r.at<double>(1, 2);
                e->R_L2R = R_L2R;
                e->T_L2R = T_L2R;
                e->Xw[0] = pt.x;
                e->Xw[1] = pt.y;
                e->Xw[2] = pt.z;

                optimizer.addEdge(e);
            }
            break;

    case 3:
            cout << "bino: optimization by lines and points" << endl;
            // add add vertex for points and edges for projection
            for (int i = 0; i < npts; i++)
            {
                const auto &pt = pts3d[i];
                const auto &l_uv = l_pts2d[i];
                const auto &r_uv = r_pts2d[i];

                g2o::EdgeSteoSE3ProjectXYZOnlyPose *e = new g2o::EdgeSteoSE3ProjectXYZOnlyPose();
                e->setVertex(0, v_se3);
                e->setMeasurement(Eigen::Vector4d(l_uv.x, l_uv.y, r_uv.x, r_uv.y));
                e->information() = Eigen::Matrix4d::Identity();
                // e->setParameterId(0, 0);
                e->setRobustKernel(new g2o::RobustKernelHuber);

                e->fx_l = _K_l.at<double>(0, 0);
                e->fy_l = _K_l.at<double>(1, 1);
                e->cx_l = _K_l.at<double>(0, 2);
                e->cy_l = _K_l.at<double>(1, 2);
                e->fx_r = _K_r.at<double>(0, 0);
                e->fy_r = _K_r.at<double>(1, 1);
                e->cx_r = _K_r.at<double>(0, 2);
                e->cy_r = _K_r.at<double>(1, 2);
                e->R_L2R = R_L2R;
                e->T_L2R = T_L2R;
                e->Xw[0] = pt.x;
                e->Xw[1] = pt.y;
                e->Xw[2] = pt.z;
                optimizer.addEdge(e);
            }

            // add vertex for lines and edges for projection
            for (int i = 0; i < nlns; i++)
            {
                const auto &ln3d = lns3d[i];
                const auto &l_ln2d = l_lns2d[i];
                const auto &r_ln2d = r_lns2d[i];

                g2o::Vector6D temp;
                Eigen::Vector2d focal_length_l(_K_l.at<double>(0, 0), _K_l.at<double>(1, 1));
                Eigen::Vector2d principle_point_l(_K_l.at<double>(0, 2), _K_l.at<double>(1, 2));
                Eigen::Vector2d focal_length_r(_K_r.at<double>(0, 0), _K_r.at<double>(1, 1));
                Eigen::Vector2d principle_point_r(_K_r.at<double>(0, 2), _K_r.at<double>(1, 2));
                temp << ln3d[0], ln3d[1], ln3d[2], ln3d[3], ln3d[4], ln3d[5];

                g2o::BinoEdgeLine *e = new g2o::BinoEdgeLine();
                e->setparams(temp, focal_length_l, focal_length_r, principle_point_l, principle_point_r);
                e->setVertex(0, v_se3);
                e->R_L2R = R_L2R;
                e->T_L2R = T_L2R;
                e->setMeasurement(g2o::Vector8d(l_ln2d[0], l_ln2d[1], l_ln2d[2], l_ln2d[3], r_ln2d[0], r_ln2d[1], r_ln2d[2], r_ln2d[3]));
                e->information() = Eigen::Matrix<double, 8, 8>::Identity();
                // e->setParameterId(0, 0);
                e->setRobustKernel(new g2o::RobustKernelHuber);
                optimizer.addEdge(e);
            }
            break;

    default:
            cerr << "wrong flags!!!  please input 1 or 2 or 3!" << endl;
            break;
    }
    // optimize
    optimizer.initializeOptimization();
    optimizer.optimize(20);
    // output result(左目相机下的位姿)
    Eigen::MatrixXd T = Eigen::Isometry3d(v_se3->estimate()).matrix();
    R = (cv::Mat_<double>(3, 3) << T(0, 0), T(0, 1), T(0, 2),
                                  T(1, 0), T(1, 1), T(1, 2),
                                  T(2, 0), T(2, 1), T(2, 2));
    oula = R2angle(R);
    // 左手坐标系定义 旋转顺序：y-x-z
    t = (cv::Mat_<double>(3, 1) << T(0, 3), T(1, 3), T(2, 3));
    t = R.t() * t;
}