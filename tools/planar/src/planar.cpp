
#define PROGRAM_NAME \
    "planar"

#define FLAGS_CASES                                                                                \
    FLAG_CASE(double, pitch, 0.0, "Pitch angle [rad]")                                             \
    FLAG_CASE(double, roll, 0.0, "Roll angle [rad]")

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
#include "se2.hpp"

void ValidateArgs() {
    RUNTIME_ASSERT(boost::filesystem::is_regular_file(ARGS_input));
}

void ValidateFlags() {
}

template<typename T>
T approximateYaw(const Eigen::Matrix<T, 3, 3>& R) {
    return std::atan2(R(1,0), R(0,0));
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
    io::Trajectory3d trajectory3d = io::read_file<io::Trajectory3d::value_type>(ARGS_input);
    RUNTIME_ASSERT(!trajectory3d.empty());

    io::Trajectory trajectory2d;

    Eigen::Isometry3d ref3d = trajectory3d.front().pose;
    io::pose_t ref2d = io::pose_t(0.0, 0.0, 0.0);

    Eigen::Matrix3d R = (Eigen::AngleAxisd(FLAGS_pitch, Eigen::Vector3d::UnitY())
                         * Eigen::AngleAxisd(FLAGS_roll, Eigen::Vector3d::UnitX())).toRotationMatrix();

    for (const io::Trajectory3d::value_type& entry3d : trajectory3d) {
        Eigen::Isometry3d pose3d = entry3d.pose;

        Eigen::Isometry3d dT = ref3d.inverse() * pose3d;
        
        Eigen::Matrix3d R_ = R * dT.linear() * R.transpose();
        Eigen::Vector3d t_ = R * dT.translation();

        io::pose_t pose2d = PoseComposition(ref2d, io::pose_t(t_(0), t_(1), approximateYaw<double>(R_)));

        io::Trajectory::value_type entry2d;
        entry2d.id = entry3d.id;
        entry2d.pose = pose2d;

        trajectory2d.push_back(entry2d);

        ref3d = pose3d;
        ref2d = pose2d;
    }

    RUNTIME_ASSERT(io::write_file<io::Trajectory::value_type>(trajectory2d, ARGS_output));

    return 0;
}
