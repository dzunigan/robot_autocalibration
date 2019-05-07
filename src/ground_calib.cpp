
#define PROGRAM_NAME \
    "ground_calib"

#define FLAGS_CASES                                                                                \
    FLAG_CASE(bool, csv, false, "Use csv input format instead of binary")                          \
    FLAG_CASE(bool, verbose, false, "Output additional information")

#define ARGS_CASES                                                                                 \
    ARG_CASE(input_file)

#include <cstdint>
#include <fstream>
#include <iostream>

// Boost
#include <boost/filesystem.hpp>

#include "estimators.hpp"
#include "refinement.hpp"

#include "util/alignment.h"
#include "util/args.hpp"
#include "util/csv.hpp"
#include "util/endian.hpp"
#include "util/macros.h"
#include "util/misc.hpp"

void ValidateArgs() {
    RUNTIME_ASSERT(boost::filesystem::is_regular_file(ARGS_input_file));
}

void ValidateFlags() {
}

inline void PrintCSV(std::ostream& out, const Eigen::Ref<const Eigen::Vector3d>& vec, int precision = 6) {
    out << std::fixed << std::setprecision(precision)
        << vec(0) << ", " << vec(1) << ", " << vec(2);
}

inline void PrintWSV(std::ostream& out, const Eigen::Ref<const Eigen::Vector3d>& vec, int precision = 6) {
    out << std::fixed << std::setprecision(precision)
        << vec(0) << " " << vec(1) << " " << vec(2);
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

    GroundCalibrationConfig config;
    if (FLAGS_csv) {
        Eigen::MatrixXd data = csv::read<double>(ARGS_input_file);
        config.Observations().reserve(data.rows());
        for (int i = 0; i < data.rows(); ++i)
            config.Observations().push_back(data.row(i).transpose());
    } else {
        std::ios::openmode mode = std::ios::in;
        if (!FLAGS_csv) mode |= std::ios::binary;

        std::ifstream input_stream(ARGS_input_file, mode);
        RUNTIME_ASSERT(input_stream.is_open());

        std::uint64_t n = ReadBinaryLittleEndian<std::uint64_t>(&input_stream);

        config.Observations().resize(n);
        for (std::uint64_t i = 0; i < n; ++i) {
            Eigen::Vector3d p;

            p(0) = ReadBinaryLittleEndian<float>(&input_stream);
            p(1) = ReadBinaryLittleEndian<float>(&input_stream);
            p(2) = ReadBinaryLittleEndian<float>(&input_stream);

            config.Observation(i) = p;
        }
    }

    std::vector<Eigen::Vector3d> models = GroundEstimator::Estimate(config.Observations());
    RUNTIME_ASSERT(!models.empty());

    config.Parameters() = models.front();

    if (FLAGS_verbose) {
        PrintHeading1("Closed-form solution");
        std::cout << "z, pitch [rad], roll [rad]" << std::endl;
        PrintCSV(std::cout, config.Parameters());
        std::cout << std::endl;
    }

    GroundCalibrationOptions options;
    options.loss_function_type = GroundCalibrationOptions::LossFunctionType::CAUCHY;
    options.loss_function_scale = 0.01;

    GroundCalibration calibration(options, config);
    calibration.Solve();

    if (FLAGS_verbose) {
        PrintHeading1("Iterative refinement");
        std::cout << "z, pitch [rad], roll [rad]" << std::endl;
        PrintCSV(std::cout, calibration.Parameters());
        std::cout << std::endl;
    } else {
        PrintWSV(std::cout, calibration.Parameters());
    }

    return 0;
}

