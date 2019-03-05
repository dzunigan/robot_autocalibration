// Calibration

// STL
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

// Boost
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/sign.hpp>

// Gflags
#include <gflags/gflags.h>

// Glog
#include <glog/logging.h>

#include "loransac.h"
#include "estimators.hpp"

#include "util/alignment.h"
#include "util/csv.hpp"
#include "util/misc.hpp"

using Calibration_RANSAC = colmap::LORANSAC<CalibrationEstimator3d, CalibrationEstimator3d>;
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(Calibration_RANSAC::Report)

#ifndef FLAGS_CASES
#define FLAGS_CASES                                                                                \
    FLAG_CASE(uint64, min_num_inliers, 6, "Minimum number of inliers to consider a link")          \
    FLAG_CASE(string, scale_ambiguous, "", "Comma separated motion ids that have scale ambiguity")
#endif

#define FLAG_CASE(type, name, val, txt) \
    DEFINE_##type(name, val, txt);

FLAGS_CASES

#undef FLAG_CASE

#ifndef PROGRAM_NAME
#define PROGRAM_NAME \
    "calibrate"
#endif

bool HelpRequired(int argc, char* argv[]) noexcept {

    const std::string help_flag("--help");
    for (int i = 1; i < argc; ++i)
        if (help_flag.compare(argv[i]) == 0) return true;
    return false;
}

void ShowHelp() noexcept {

    std::cerr << "Usage: " << PROGRAM_NAME << " [options]";
    std::cerr << " <reference> <sensor> ";
    std::cerr << std::endl;

    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;

#define FLAG_CASE(type, name, val, txt)                                                            \
    std::cerr << "  --" << #name << ": " << txt << std::endl;                                      \
    std::cerr << "        " << "(type: " << #type << ", default: " << #val << ")" << std::endl;

    FLAGS_CASES

#undef FLAG_CASE

    std::cerr << "  --help: Displays this message" << std::endl;
    std::cerr << std::endl;
}

inline void PrintCSV(std::ostream& out, const Eigen::Ref<const Eigen::Vector6d>& vec, int precision = 6) {
    out << std::fixed << std::setprecision(precision)
        << vec(0) << ", " << vec(1) << ", " << vec(2) << ", " << vec(3) << ", " << vec(4) << ", " << vec(5);
}

void ValidateFlags() {
    if (!FLAGS_scale_ambiguous.empty()) {
        CHECK(!VectorContainsValue(CSVToVector<int>(FLAGS_scale_ambiguous), 1));
        CHECK(!VectorContainsDuplicateValues(CSVToVector<int>(FLAGS_scale_ambiguous)));
    }
}

void ValidateArgs(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i)
        CHECK(boost::filesystem::is_regular_file(argv[i])) << "Invalid input file:" << std::endl
                                                           << argv[i];
}

int main(int argc, char* argv[]) {

    // Check input args
    // -------------------

    // Handle help flag
    if (HelpRequired(argc, argv)) {
        ShowHelp();
        return 0;
    }

    // Initialize Google's logging library
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = true;

    // Parse input flags
    gflags::ParseCommandLineNonHelpFlags(&argc, &argv, true);

    // Two sensors required
    if (argc != 3) {
        ShowHelp();
        return -1;
    }

    // Check input arguments
    ValidateFlags();
    ValidateArgs(argc, argv);

    // Load sensor observations
    // -------------------
    std::vector<Eigen::Vector6d> reference;
    {
        Eigen::MatrixXd data = csv::read<double>(argv[1], ' ');
        CHECK_EQ(data.cols(), 7);

        for (int k = 0; k < data.rows(); ++k) {
            Eigen::Vector3d t_vec = data.block<1, 3>(k, 0).transpose();
            Eigen::Vector4d q_vec = data.block<1, 4>(k, 4).transpose();

            Eigen::Quaterniond q(q_vec(3), q_vec(0), q_vec(1), q_vec(2));
            q.normalize();

            Eigen::AngleAxisd angle_axis(q);

            Eigen::Vector6d v;
            v.head<3>() = t_vec;
            v.tail<3>() = angle_axis.angle() * angle_axis.axis();

            reference.push_back(v);
        }
    }

    std::vector<Eigen::Vector6d> sensor;
    {
        Eigen::MatrixXd data = csv::read<double>(argv[2], ' ');
        CHECK_EQ(data.cols(), 7);

        for (int k = 0; k < data.rows(); ++k) {
            Eigen::Vector3d t_vec = data.block<1, 3>(k, 0).transpose();
            Eigen::Vector4d q_vec = data.block<1, 4>(k, 4).transpose();

            Eigen::Quaterniond q(q_vec(3), q_vec(0), q_vec(1), q_vec(2));
            q.normalize();

            Eigen::AngleAxisd angle_axis(q);

            Eigen::Vector6d v;
            v.head<3>() = t_vec;
            v.tail<3>() = angle_axis.angle() * angle_axis.axis();

            sensor.push_back(v);
        }
    }

    // All sensor must have the same number of observations
    CHECK_EQ(reference.size(), sensor.size());

    // Auto-initialization
    // -------------------

    colmap::RANSACOptions ransac_options;
    ransac_options.max_error = 0.1;
    ransac_options.min_inlier_ratio = 0.25;
    ransac_options.confidence = 0.99999;
    ransac_options.min_num_trials = 100;
    ransac_options.max_num_trials = 1000;

    Calibration_RANSAC::Report ransac_report;

    Calibration_RANSAC calibration_ransac(ransac_options);
    ransac_report = calibration_ransac.Estimate(reference, sensor);

    // Check inlier set
    CHECK(ransac_report.inlier_mask.size() > ransac_options.min_inlier_ratio * sensor.size()) << "Not enough inliers!";

    PrintCSV(std::cout, ransac_report.model);
    std::cout << std::endl;

    return 0;
}
