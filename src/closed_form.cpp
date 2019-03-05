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

#include "refinement.hpp"

#include "util/alignment.h"
#include "util/csv.hpp"
#include "util/misc.hpp"

using Calibration_RANSAC = colmap::LORANSAC<CalibrationEstimator, CalibrationEstimator>;
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(Calibration_RANSAC::Report)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(Sensor)

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
    std::cerr << " <motions1> <motions2> ...";
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

inline void PrintCSV(std::ostream& out, const Eigen::Ref<const Eigen::Vector4d>& vec, int precision = 6) {
    out << std::fixed << std::setprecision(precision)
        << vec(0) << ", " << vec(1) << ", " << vec(2) << ", " << vec(3);
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

    // Check at runtime that size_t can contain all sensor ids without information loss
    CHECK_EQ(static_cast<size_t>(kMaxNumSensors), kMaxNumSensors);

    // Parse input flags
    gflags::ParseCommandLineNonHelpFlags(&argc, &argv, true);

    // At least two sensors required
    if (argc < 3) {
        ShowHelp();
        return -1;
    }

    // Check input arguments
    ValidateFlags();
    ValidateArgs(argc, argv);

    // Load sensor observations
    // -------------------
    std::vector<Sensor> sensors(argc-1);
    for (int i = 1; i < argc; ++i) {
        Eigen::MatrixXd data = csv::read<double>(argv[i], ' ');

        Sensor& sensor = sensors[i-1];
        sensor.SetSensorId(i-1);

        for (int k = 0; k < data.rows(); ++k)
            sensor.Observations().push_back(data.row(k));
    }

    // Set reference sensor
    const Sensor& reference = sensors.front();
    const size_t N = reference.Observations().size();

    // All sensor must have the same number of observations
    for (const Sensor& sensor : sensors)
        CHECK_EQ(sensor.Observations().size(), N);

    // Auto-initialization
    // -------------------
    ObservationDatabase observations; // RANSAC inliers

    colmap::RANSACOptions ransac_options;
    ransac_options.max_error = 0.03;
    ransac_options.min_inlier_ratio = 0.15;
    ransac_options.confidence = 0.99999;
    ransac_options.min_num_trials = 1000;
    ransac_options.max_num_trials = 5000;

    std::vector<Calibration_RANSAC::Report> reports(sensors.size() - 1);
    for (size_t i = 1; i < sensors.size(); ++i) {
        Sensor& sensor = sensors[i];
        Calibration_RANSAC::Report& ransac_report = reports[i-1];

        Calibration_RANSAC calibration_ransac(ransac_options);
        ransac_report = calibration_ransac.Estimate(reference.Observations(), sensor.Observations());

        CHECK_EQ(sensor.Observations().size(), ransac_report.inlier_mask.size());
        for (size_t k = 0; k < sensor.Observations().size(); ++k) {
            if (ransac_report.inlier_mask.at(k))
                observations.AddIndexPair(reference.SensorId(), sensor.SensorId(), k, k);
        }

        // Check inlier set
        CHECK(ransac_report.inlier_mask.size() > ransac_options.min_inlier_ratio * sensor.Observations().size()) << "Not enough inliers!";
        //sensor.Calibration() = ransac_report.model;
        PrintCSV(std::cout, ransac_report.model);
        std::cout << std::endl;
    }

    return 0;
}
