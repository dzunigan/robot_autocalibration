
#define PROGRAM_NAME \
    "scale_trajectory"

#define FLAGS_CASES                                                                                \
    FLAG_CASE(double, scale, 1.0, "Scale factor")

#define ARGS_CASES                                                                                 \
    ARG_CASE(input)                                                                                \
    ARG_CASE(output)

// STL
#include <iostream>

// Boost
#include <boost/filesystem.hpp>

// Eigen
#include <Eigen/Geometry>

#include "args.hpp"
#include "macros.h"
#include "io3d.hpp"

void ValidateArgs() {
    RUNTIME_ASSERT(boost::filesystem::is_regular_file(ARGS_input));
}

void ValidateFlags() {
    RUNTIME_ASSERT(FLAGS_scale > 0.0);
}

int main(int argc, char* argv[]) {

    // Handle help flag
    if (args::HelpRequired(argc, argv)) {
        args::ShowHelp();
        return 0;
    }

    // Parse input flags
    gflags::ParseCommandLineNonHelpFlags(&argc, &argv, true);

    // Check number of args
    if (argc-1 != args::NumArgs()) {
        args::ShowHelp();
        return -1;
    }

    // Parse input args
    args::ParseCommandLineArgs(argc, argv);

    // Check input arguments
    ValidateFlags();
    ValidateArgs();

    // Read input trajectories
    io::Trajectory3d trajectory = io::read_file<io::Trajectory3d::value_type>(ARGS_input);
    RUNTIME_ASSERT(!trajectory.empty());

    io::Trajectory3d trajectory_scaled;
    for (const io::Trajectory3d::value_type& entry : trajectory) {
        Eigen::Isometry3d pose = entry.pose;
        pose.translation() *= FLAGS_scale;

        io::Trajectory3d::value_type entry_scaled;
        entry_scaled.id = entry.id;
        entry_scaled.pose = pose;

        trajectory_scaled.push_back(entry_scaled);
    }

    RUNTIME_ASSERT(io::write_file<io::Trajectory3d::value_type>(trajectory_scaled, ARGS_output));

    return 0;
}
