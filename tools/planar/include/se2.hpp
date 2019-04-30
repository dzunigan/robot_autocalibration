
#pragma once

#include <cmath>

// Eigen
#include <Eigen/Core>

#include "io.hpp"

// Normalizes the angle in radians between [-pi and pi).
inline double NormalizeAngle(const double angle_radians) {
    const double two_pi(2.0 * EIGEN_PI);
    return angle_radians -
            two_pi * std::floor((angle_radians + EIGEN_PI) / two_pi);
}

inline io::pose_t PoseComposition(const io::pose_t& a, const io::pose_t& b) {

    const double cos_theta = std::cos(a.theta);
    const double sin_theta = std::sin(a.theta);

    io::pose_t c;
    c.tx = a.tx + b.tx * cos_theta - b.ty * sin_theta;
    c.ty = a.ty + b.tx * sin_theta + b.ty * cos_theta;
    c.theta = NormalizeAngle(a.theta + b.theta);

    return c;
}

inline io::pose_t PoseInverse(const io::pose_t& a) {

    const double cos_theta = std::cos(a.theta);
    const double sin_theta = std::sin(a.theta);

    io::pose_t inv;
    inv.tx = -a.tx * cos_theta - a.ty * sin_theta;
    inv.ty =  a.tx * sin_theta - a.ty * cos_theta;
    inv.theta = NormalizeAngle(-a.theta);

    return inv;
}
