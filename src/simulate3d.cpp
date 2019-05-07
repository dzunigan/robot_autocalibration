
#define PROGRAM_NAME \
    "simulate3d"

#define FLAGS_CASES                                                                                \
    FLAG_CASE(double, lin_std, 0.005, "Standard deviation for translation [m]")                     \
    FLAG_CASE(double, ang_std, 0.03, "Standard deviation for rotation [rad]")

#define ARGS_CASES                                                                                 \
    ARG_CASE(output1)                                                                              \
    ARG_CASE(output2)

#include <iostream>

#include "util/args.hpp"
#include "util/macros.h"
#include "util/io3d.hpp"
#include "util/random.h"

#include "sim2.hpp"

#include <Eigen/Geometry>

void ValidateArgs() {
}

void ValidateFlags() {
    RUNTIME_ASSERT(FLAGS_lin_std >= 0.0);
    RUNTIME_ASSERT(FLAGS_ang_std >= 0.0);
}

Eigen::Isometry3d simulate_error(const Eigen::Isometry3d& p) {
    using namespace colmap;

    Eigen::Isometry3d e = Eigen::Isometry3d::Identity();

    e.linear() = (Eigen::AngleAxisd(RandomGaussian(0.0, FLAGS_ang_std) * EIGEN_PI, Eigen::Vector3d::UnitZ()) *
                  Eigen::AngleAxisd(RandomGaussian(0.0, FLAGS_ang_std) * EIGEN_PI, Eigen::Vector3d::UnitY()) *
                  Eigen::AngleAxisd(RandomGaussian(0.0, FLAGS_ang_std) * EIGEN_PI, Eigen::Vector3d::UnitX()))
                  .toRotationMatrix();
 
    e.translation() = Eigen::Vector3d(RandomGaussian(0.0, FLAGS_lin_std),
                                      RandomGaussian(0.0, FLAGS_lin_std),
                                      RandomGaussian(0.0, FLAGS_lin_std));

    return p * e;
}

Eigen::Isometry3d simulate_odometry_error(const Eigen::Isometry3d& p) {
    using namespace colmap;

    Eigen::Isometry3d e = Eigen::Isometry3d::Identity();

    e.linear() = Eigen::AngleAxisd(RandomGaussian(0.0, FLAGS_ang_std) * EIGEN_PI, Eigen::Vector3d::UnitZ())
                 .toRotationMatrix();
 
    e.translation() = Eigen::Vector3d(RandomGaussian(0.0, FLAGS_lin_std),
                                      RandomGaussian(0.0, FLAGS_lin_std),
                                      0.0);

    return p * e;
}

void add_pose(io::Trajectory3d& traj, io::timestamp_t t, const Eigen::Isometry3d& p) {
    io::trajectory3d_t<io::timestamp_t> entry;
    entry.id = t;
    entry.pose = io::pose3d_t(p);

    traj.push_back(entry);
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

    const double linear_vel = 0.5; // m/s
    const double angular_vel = 0.5; // rad/s

    const double scale_factor = 0.5;
    Eigen::Isometry3d rel;
    rel.linear() = (Eigen::AngleAxisd(- 1.0 / 2.0 * EIGEN_PI, Eigen::Vector3d::UnitZ())
                    * Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitY())
                    * Eigen::AngleAxisd(-3.0 / 4.0 * EIGEN_PI, Eigen::Vector3d::UnitX())).toRotationMatrix();
    rel.translation() = Eigen::Vector3d(0.1, 0.5, 1.0);

    Eigen::Isometry3d last_pose1;
    last_pose1.linear() = Eigen::Matrix3d::Identity();
    last_pose1.translation() = Eigen::Vector3d::Zero();

    Eigen::Isometry3d last_pose2;
    last_pose2.linear() = Eigen::Matrix3d::Identity();
    last_pose2.translation() = Eigen::Vector3d::Zero();

    Eigen::Isometry3d dT;
    dT.linear() = Eigen::AngleAxisd(angular_vel, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    dT.translation() = Eigen::Vector3d(linear_vel, 0.0, 0.0);

    io::timestamp_t t = 0;
    io::Trajectory3d t1, t2;

    for (int k = 0; k < 3; ++k) {
      for (int i = 0; i < 12; ++i) {
          last_pose1 = simulate_odometry_error(last_pose1 * dT);
          last_pose2 = simulate_odometry_error(last_pose2 * (rel.inverse() * dT * rel));

          t += 1;

          last_pose2.translation() *= scale_factor;

          add_pose(t1, t, last_pose1);
          add_pose(t2, t, last_pose2);

          last_pose2.translation() /= scale_factor;
      }

      dT.linear() = Eigen::AngleAxisd(-angular_vel, Eigen::Vector3d::UnitZ()).toRotationMatrix();

      for (int i = 0; i < 12; ++i) {
          last_pose1 = simulate_odometry_error(last_pose1 * dT);
          last_pose2 = simulate_odometry_error(last_pose2 * (rel.inverse() * dT * rel));

          t += 1;

          last_pose2.translation() *= scale_factor;

          add_pose(t1, t, last_pose1);
          add_pose(t2, t, last_pose2);

          last_pose2.translation() /= scale_factor;
      }

      last_pose1 = simulate_odometry_error(last_pose1 * dT);
      last_pose2 = simulate_odometry_error(last_pose2 * (rel.inverse() * dT * rel));

      t += 1;

      last_pose2.translation() *= scale_factor;

      add_pose(t1, t, last_pose1);
      add_pose(t2, t, last_pose2);

      last_pose2.translation() /= scale_factor;

      dT.linear() = Eigen::AngleAxisd(angular_vel, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    }

    // Save trajectories
    RUNTIME_ASSERT(io::write_file<io::Trajectory3d::value_type>(t1, ARGS_output1));
    RUNTIME_ASSERT(io::write_file<io::Trajectory3d::value_type>(t2, ARGS_output2));

    return 0;
}
