
#define PROGRAM_NAME \
    "from2dto3d"

#define FLAGS_CASES                                                                                \

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
#include "io.hpp"
#include "io3d.hpp"

void ValidateArgs() {
    RUNTIME_ASSERT(boost::filesystem::is_regular_file(ARGS_input));
}

void ValidateFlags() {
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

    // Read input trajectory
    io::Trajectory trajectory2d = io::read_file<io::Trajectory::value_type>(ARGS_input);
    RUNTIME_ASSERT(!trajectory2d.empty());

    io::Trajectory3d trajectory3d;
    for (const io::Trajectory::value_type& entry2d : trajectory2d) {
        Eigen::Isometry3d T;
        T.linear() = Eigen::AngleAxisd(entry2d.pose.theta, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        T.translation() = Eigen::Vector3d(entry2d.pose.tx, entry2d.pose.ty, 0.0);

        io::pose3d_t pose3d = io::pose3d_t(T);

        io::Trajectory3d::value_type entry3d;
        entry3d.id = entry2d.id;
        entry3d.pose = pose3d;

        trajectory3d.push_back(entry3d);
    }

    RUNTIME_ASSERT(io::write_file<io::Trajectory3d::value_type>(trajectory3d, ARGS_output));

    return 0;
}
