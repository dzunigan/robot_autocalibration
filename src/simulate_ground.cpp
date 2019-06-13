
#define PROGRAM_NAME \
    "simulate_ground"

#define FLAGS_CASES                                                                                \
    FLAG_CASE(double, point_std, 0.01, "Standard deviation for 3D points [m]")

#define ARGS_CASES                                                                                 \
    ARG_CASE(output)

#include <iostream>

#include "util/args.hpp"
#include "util/csv.hpp"
#include "util/macros.h"
#include "util/math.h"
#include "util/io3d.hpp"
#include "util/random.h"

#include "sim2.hpp"

#include <Eigen/Geometry>

void ValidateArgs() {
}

void ValidateFlags() {
    RUNTIME_ASSERT(FLAGS_point_std >= 0.0);
}

Eigen::Vector3d simulate_error(const Eigen::Vector3d& p) {
    using namespace colmap;
    Eigen::Vector3d e(0.0, 0.0, RandomGaussian(0.0, FLAGS_point_std));

    return p + e;
}

int main(int argc, char* argv[]) {

    // Handle help flag
    if (args::HelpRequired(argc, argv)) {
        args::ShowHelp();
        return 0;
    }

    // Parse input flags
    args::ParseCommandLineNonHelpFlags(&argc, &argv, true);

    // Check number of args
    if (argc-1 != args::NumArgs()) {
        args::ShowHelp();
        return -1;
    }

    // Parse input args
    args::ParseCommandLineArgs(argc, argv);

    // Validate input arguments
    ValidateArgs();
    ValidateFlags();

    const double scale_factor = 0.5;
    Eigen::Isometry3d rel;
    rel.linear() = (Eigen::AngleAxisd(- 1.0 / 2.0 * EIGEN_PI, Eigen::Vector3d::UnitZ())
                    * Eigen::AngleAxisd(1.0 / 12.0, Eigen::Vector3d::UnitY())
                    * Eigen::AngleAxisd(-3.0 / 4.0 * EIGEN_PI, Eigen::Vector3d::UnitX())).toRotationMatrix();
    rel.translation() = Eigen::Vector3d(0.5, 0.1, 1.0);

    Eigen::Vector4d plane(0.0, 0.0, 1.0, 0.0);
    plane = hesseNormalForm(plane);
    plane = transform_plane(plane, rel.inverse());

    // Camera parameters
    const unsigned int width = 320;
    const unsigned int height = 240;
    const double f = 285.2;

    const double cx = width / 2 - 0.5;
    const double cy = height / 2 - 0.5;
    const double inv_fx = 1.0 / f;
    const double inv_fy = 1.0 / f;

    Eigen::Isometry3d ref;
    ref.linear() = Eigen::Matrix3d::Identity();
    ref.translation() = Eigen::Vector3d::Zero();

    Eigen::MatrixXd points(width * height, 3);

    int idx = 0;
    for (unsigned int i = 0; i < height; ++i) {
        for (unsigned int j = 0; j < width; ++j) {
            const double coef_x = (j - cx) * inv_fx;
            const double coef_y = (i - cy) * inv_fy;

            Eigen::Vector3d l(coef_x, coef_y, 1.0);

            const double z = -plane(3) / plane.head<3>().dot(l);

            Eigen::Vector3d p = (l * z);

            points.row(idx) = simulate_error(p).transpose();
            points.row(idx) *= scale_factor;
            idx++;
        }
    }

    RUNTIME_ASSERT(csv::write(points, ARGS_output));

    return 0;
}
