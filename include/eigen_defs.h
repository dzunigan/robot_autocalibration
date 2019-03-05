#ifndef CALIBRATION_EIGEN_DEFS_H_
#define CALIBRATION_EIGEN_DEFS_H_

// STL
#include <vector>

// Eigen
#include <Eigen/Core>
#include <Eigen/StdVector>

#include "util/alignment.h"

namespace Eigen {
    template<typename T>
    using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;

    using Matrix5d = Eigen::Matrix<double, 5, 5>;
    using Vector5d = Eigen::Matrix<double, 5, 1>;

    using Matrix6d = Eigen::Matrix<double, 6, 6>;
    using Vector6d = Eigen::Matrix<double, 6, 1>;
}

#endif // CALIBRATION_EIGEN_DEFS_H_
