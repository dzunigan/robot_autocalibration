
#pragma once

// Cere
#include <ceres/jet.h>

// Eigen
#include <Eigen/Core>

#include "normalize_angle.h"

template <typename T>
inline Eigen::Matrix<T, 2, 2> RotationMatrix2D(T yaw_radians) {
    const T cos_yaw = ceres::cos(yaw_radians);
    const T sin_yaw = ceres::sin(yaw_radians);

    Eigen::Matrix<T, 2, 2> rotation;
    rotation << cos_yaw, -sin_yaw, sin_yaw, cos_yaw;
    return rotation;
}

template<typename T>
inline Eigen::Matrix<T, 4, 1> PoseComposition(const Eigen::Ref<const Eigen::Matrix<T, 4, 1>> &a, const Eigen::Ref<const Eigen::Matrix<T, 4, 1>> &b) {

    const T cos_theta = ceres::cos(a(2));
    const T sin_theta = ceres::sin(a(2));

    Eigen::Matrix<T, 4, 1> c;
    //c.head<2>() = (T(1.0) / b(3)) * a.head<2>() + RotationMatrix2D(a(2)) * b.head<2>();

    c(0) = (T(1.0) / b(3)) * a(0) + b(0) * cos_theta - b(1) * sin_theta;
    c(1) = (T(1.0) / b(3)) * a(1) + b(0) * sin_theta + b(1) * cos_theta;
    c(2) = NormalizeAngle(a(2) + b(2));
    c(3) = a(3) * b(3);

    return c;
}

template<typename T>
inline Eigen::Matrix<T, 4, 1> PoseInverse(const Eigen::Ref<const Eigen::Matrix<T, 4, 1>> &a) {

    const T cos_theta = ceres::cos(a(2));
    const T sin_theta = ceres::sin(a(2));

    Eigen::Matrix<T, 4, 1> inv;
    //inv.head(2) = a(3) * RotationMatrix2D(a(2)).transpose() * a.head(2);

    inv(0) = -a(3) * (a(0) * cos_theta + a(1) * sin_theta);
    inv(1) = -a(3) * (-a(0) * sin_theta + a(1) * cos_theta);
    inv(2) = -a(2);
    inv(3) = T(1.0) / a(3);

    return inv;
}
