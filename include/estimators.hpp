
// Eigen
#include <Eigen/Core>

#include "eigen_defs.h"

/*
class RotationAndScaleEstimator {
public:
    typedef Eigen::Vector3d X_t;
    typedef Eigen::Vector3d Y_t;
    typedef Eigen::Vector2d M_t; // (scale factor, angle)

    // The minimum number of samples needed to estimate a model.
    static const int kMinNumSamples = 1;

    // Estimate extrinsic calibration solutions from a set of
    // corresponding incremental motions.
    //
    // At least a single correspondence needed
    //
    // @param s1  First set of corresponding incremental motions.
    // @param s2  Second set of corresponding incremental motions.
    //
    // @return    Relative pose from s1 to s2
    static std::vector<M_t> Estimate(const std::vector<X_t>& s1,
                                     const std::vector<Y_t>& s2);

    // Calculate the residuals of a set of corresponding incremental motions
    // and a given extrinsic calibration.
    //
    // Residuals are defined as the squared Sampson error.
    //
    // @param s1         First set of corresponding incremental motions.
    // @param s2         Second set of corresponding incremental motions.
    // @param c          Extrinsic sensor calibration.
    // @param residuals  Output vector of residuals.
    static void Residuals(const std::vector<X_t>& s1,
                        const std::vector<Y_t>& s2, const M_t& c,
                        std::vector<double>* residuals);
};
*/

class CalibrationEstimator {
public:
    typedef Eigen::Vector3d X_t;
    typedef Eigen::Vector3d Y_t;
    typedef Eigen::Vector4d M_t;

    // The minimum number of samples needed to estimate a model.
    static const int kMinNumSamples = 2;

    // Estimate extrinsic calibration solutions from a set of
    // corresponding incremental motions.
    //
    // At least two correspondences are needed
    //
    // @param x   First set of incremental motions (metrically accurate, e.g. wheel odometry).
    // @param y   Second set of incremental motions (possibly scale ambiguous, e.g. monocular vo).
    //
    // @return    Relative pose of sensor observing y's trajectory to x's sensor
    static std::vector<M_t> Estimate(const std::vector<X_t>& x,
                                     const std::vector<Y_t>& y);

    // Calculate the residuals of a set of corresponding incremental motions
    // given an extrinsic calibration model.
    //
    // Residuals are defined as the squared translational error in the x's metric space.
    //
    // @param x          First set of incremental motions (metrically accurate, e.g. wheel odometry).
    // @param y          Second set of incremental motions (possibly scale ambiguous, e.g. monocular vo).
    // @param m          Model representing the extrinsic calibration between x's and y's sensors.
    // @param residuals  Output vector of residuals.
    static void Residuals(const std::vector<X_t>& x,
                          const std::vector<Y_t>& y, const M_t& m,
                          std::vector<double>* residuals);
};

class CalibrationEstimator2 {
public:
    typedef Eigen::Vector3d X_t;
    typedef Eigen::Vector3d Y_t;
    typedef Eigen::Vector4d M_t;

    // The minimum number of samples needed to estimate a model.
    static const int kMinNumSamples = 2;

    // Estimate extrinsic calibration solutions from a set of
    // corresponding incremental motions.
    //
    // At least two correspondences are needed
    //
    // @param x   First set of incremental motions (metrically accurate, e.g. wheel odometry).
    // @param y   Second set of incremental motions (possibly scale ambiguous, e.g. monocular vo).
    //
    // @return    Relative pose of sensor observing y's trajectory to x's sensor
    static std::vector<M_t> Estimate(const std::vector<X_t>& x,
                                     const std::vector<Y_t>& y);

    // Calculate the residuals of a set of corresponding incremental motions
    // given an extrinsic calibration model.
    //
    // Residuals are defined as the squared translational error in the x's metric space.
    //
    // @param x          First set of incremental motions (metrically accurate, e.g. wheel odometry).
    // @param y          Second set of incremental motions (possibly scale ambiguous, e.g. monocular vo).
    // @param m          Model representing the extrinsic calibration between x's and y's sensors.
    // @param residuals  Output vector of residuals.
    static void Residuals(const std::vector<X_t>& x,
                          const std::vector<Y_t>& y, const M_t& m,
                          std::vector<double>* residuals);
};

class CalibrationEstimator3d {
public:
    typedef Eigen::Vector6d X_t; // x, y, z, ax, ay, az (axis angle)
    typedef Eigen::Vector6d Y_t; // x, y, z, ax, ay, az (axis angle)
    typedef Eigen::Vector6d M_t; // x, y, yaw, pitch, roll, scale

    // The minimum number of samples needed to estimate a model.
    static const int kMinNumSamples = 2;

    // Estimate extrinsic calibration solutions from a set of
    // corresponding incremental motions.
    //
    // At least two correspondences are needed
    //
    // @param x   First set of incremental motions (metrically accurate, e.g. wheel odometry).
    // @param y   Second set of incremental motions (possibly scale ambiguous, e.g. monocular vo).
    //
    // @return    Relative pose of sensor observing y's trajectory to x's sensor
    static std::vector<M_t> Estimate(const std::vector<X_t>& x,
                                     const std::vector<Y_t>& y);

    // Calculate the residuals of a set of corresponding incremental motions
    // given an extrinsic calibration model.
    //
    // Residuals are defined as the squared translational error in the x's metric space.
    //
    // @param x          First set of incremental motions (metrically accurate, e.g. wheel odometry).
    // @param y          Second set of incremental motions (possibly scale ambiguous, e.g. monocular vo).
    // @param m          Model representing the extrinsic calibration between x's and y's sensors.
    // @param residuals  Output vector of residuals.
    static void Residuals(const std::vector<X_t>& x,
                          const std::vector<Y_t>& y, const M_t& m,
                          std::vector<double>* residuals);

private:
    static bool estimateRyx(const std::vector<X_t>& x,
                                     const std::vector<Y_t>& y, Eigen::Matrix3d& R_yx);

};