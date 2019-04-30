
#define PROGRAM_NAME \
    "simulate"

#define FLAGS_CASES                                                                                \
    FLAG_CASE(double, lin_std, 0.01, "Standard deviation for translation [m]")                     \
    FLAG_CASE(double, ang_std, 0.03, "Standard deviation for rotation [rad]")

#define ARGS_CASES                                                                                 \
    ARG_CASE(output1)                                                                              \
    ARG_CASE(output2)

#include <iostream>

#include "util/args.hpp"
#include "util/macros.h"
#include "util/io.hpp"
#include "util/random.h"

#include "sim2.hpp"

void ValidateArgs() {
}

void ValidateFlags() {
    RUNTIME_ASSERT(FLAGS_lin_std >= 0.0);
    RUNTIME_ASSERT(FLAGS_ang_std >= 0.0);
}

Eigen::Vector4d simulate_error(const Eigen::Vector4d& p) {
    using namespace colmap;
    Eigen::Vector4d e(RandomGaussian(0.0, FLAGS_lin_std), RandomGaussian(0.0, FLAGS_lin_std), RandomGaussian(0.0, FLAGS_ang_std), 1.0);

    return PoseComposition<double>(p, e);
}

void add_pose(io::Trajectory& traj, io::timestamp_t t, const Eigen::Vector4d& p) {
    io::trajectory_t<io::timestamp_t> entry;
    entry.id = t;
    entry.pose = io::pose_t(p(0), p(1), p(2));

    traj.push_back(entry);
}

int main(int argc, char* argv[]) {

    colmap::SetPRNGSeed(1023);

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

    const double linear_vel = 0.5; // m/s
    const double angular_vel = 0.5; // rad/s

    Eigen::Vector4d rel(0.5, 0.1, -0.2, 1);

    Eigen::Vector4d ref(0, 0, 0, 1);
    Eigen::Vector4d dT;

    dT = Eigen::Vector4d(linear_vel, 0, angular_vel, 1);

    io::timestamp_t t = 0;
    io::Trajectory t1, t2;

    add_pose(t1, t, simulate_error(ref));
    add_pose(t2, t, simulate_error(PoseComposition<double>(ref, rel)));

    for (int i = 0; i < 12; ++i) {
        Eigen::Vector4d pose;
        pose = PoseComposition<double>(ref, dT);

        t += 1;

        add_pose(t1, t, simulate_error(pose));
        add_pose(t2, t, simulate_error(PoseComposition<double>(pose, rel)));

        ref = pose;
    }

    dT = Eigen::Vector4d(linear_vel, 0, -angular_vel, 1);

    for (int i = 0; i < 12; ++i) {
        Eigen::Vector4d pose;
        pose = PoseComposition<double>(ref, dT);

        t += 1;

        add_pose(t1, t, simulate_error(pose));
        add_pose(t2, t, simulate_error(PoseComposition<double>(pose, rel)));

        ref = pose;
    }

    ref = PoseComposition<double>(ref, dT);

    t += 1;

    add_pose(t1, t, simulate_error(ref));
    add_pose(t2, t, simulate_error(PoseComposition<double>(ref, rel)));

    RUNTIME_ASSERT(io::write_file<io::Trajectory::value_type>(t1, ARGS_output1));
    RUNTIME_ASSERT(io::write_file<io::Trajectory::value_type>(t2, ARGS_output2));

    return 0;
}
