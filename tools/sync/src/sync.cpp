
// STL
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <utility>

// Boost
#include <boost/filesystem.hpp>

// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// Glog
#include <glog/logging.h>

#include "io.hpp"
#include "se2.hpp"

#ifndef FLAGS_CASES
#define FLAGS_CASES                                                                                \
    FLAG_CASE(uint64, s, 0, "Number of frames to skip in reference trajectory")                    \
    FLAG_CASE(string, output_dir, "", "Output directory")
#endif

#define FLAG_CASE(type, name, val, txt) \
    DEFINE_##type(name, val, txt);

FLAGS_CASES

#undef FLAG_CASE

#ifndef PROGRAM_NAME
#define PROGRAM_NAME \
    "sync"
#endif

bool HelpRequired(int argc, char* argv[]) noexcept {

    const std::string help_flag("--help");
    for (int i = 1; i < argc; ++i)
        if (help_flag.compare(argv[i]) == 0) return true;
    return false;
}

void ShowHelp() noexcept {

    std::cerr << "Usage: " << PROGRAM_NAME << " [options]";
    std::cerr << " <trajectory1> <trajectory2> ...";
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

void ValidateFlags() {
    if (!FLAGS_output_dir.empty())
        CHECK(boost::filesystem::is_directory(FLAGS_output_dir));
}

void ValidateArgs(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i)
        CHECK(boost::filesystem::is_regular_file(argv[i])) << "Invalid trayectory file:" << std::endl
                                                           << argv[i];
}

io::pose_t relative_pose(io::pose_t a, io::pose_t b) {
    io::pose_t relative;

    relative = PoseComposition(PoseInverse(a), b);

    return relative;
}

io::pose_t interpolate_pose(io::pose_t a, io::pose_t b, double t) {
    io::pose_t interpolated;

    // Linear interpolation
    interpolated.tx = a.tx + t*(b.tx - a.tx);
    interpolated.ty = a.ty + t*(b.ty - a.ty);
    interpolated.theta = NormalizeAngle(a.theta + t*NormalizeAngle(b.theta - a.theta));

    return interpolated;
}

inline bool cmp(const io::trajectory_t<double>& lhs, const double rhs) {
    return (lhs.id < rhs);
}

io::pose_t sample(const io::Trajectory trajectory, double timestamp) {
    io::Trajectory::const_iterator cit_b = std::lower_bound(trajectory.cbegin(), trajectory.cend(), timestamp, cmp);
    CHECK(cit_b != trajectory.cend());

    if (cit_b == trajectory.cbegin()) {
        CHECK_EQ(timestamp, cit_b->id);

        return cit_b->pose;
    } else {
        io::Trajectory::const_iterator cit_a = std::prev(cit_b);

        // Interpolate
        const double t = (timestamp - cit_a->id) / (cit_b->id - cit_a->id);
        return interpolate_pose(cit_a->pose, cit_b->pose, t);
    }
}

int main(int argc, char* argv[]) {

    // Handle help flag
    if (HelpRequired(argc, argv)) {
        ShowHelp();
        return 0;
    }

    // Initialize Glog library
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = true;

    // Parse input flags
    gflags::ParseCommandLineNonHelpFlags(&argc, &argv, true);

    if (argc < 3) {
        ShowHelp();
        return -1;
    }

    // Check input arguments
    ValidateFlags();
    ValidateArgs(argc, argv);

    // Read input trajectories
    std::vector<double> start_times(argc-1);
    std::vector<io::Trajectory> trajectories(argc-1);
    for (int i = 1; i < argc; ++i) {
        io::Trajectory trajectory = io::read_file<io::trajectory_t<double>>(argv[i]);
        CHECK(!trajectory.empty()) << "Empty trajectory file: " << argv[i];

        start_times[i-1] = trajectory.front().id;
        trajectories[i-1] = std::move(trajectory);
    }

    const io::Trajectory& ref = trajectories.front();

    // Initialize
    using Poses = std::vector<io::pose_t>;
    Poses last_poses(argc-1);
    size_t k = 0;
    for (; k < ref.size(); ++k) {
        double timestamp = ref[k].id;
        last_poses[0] = ref[k].pose;
        bool f = true;

        for (int i = 2; i < argc; ++i) {
            const io::Trajectory& trajectory = trajectories[i-1];

            // find lb
            io::Trajectory::const_iterator cit_b = std::lower_bound(trajectory.cbegin(), trajectory.cend(), timestamp, cmp);
            if ((cit_b == trajectory.begin() && timestamp != cit_b->id) ||
                cit_b == trajectory.cend()) {
                f = false;
                break;
            }

            // first_poses
            last_poses[i-1] = sample(trajectory, timestamp);
        }

        if (f) break;
    }

    CHECK(k < ref.size());

    bool has_next;
    std::vector<Poses> incremental_poses(argc-1);
    do {
        size_t next_k = k + FLAGS_s + 1;
        if (next_k >= ref.size()) break;

        Poses p;
        {
            io::pose_t incremental = relative_pose(last_poses[0], ref[next_k].pose);
            p.push_back(incremental);
            last_poses[0] = ref[next_k].pose;
        }

        double timestamp = ref[next_k].id;
        for (int i = 2; i < argc; ++i) {
            const io::Trajectory& trajectory = trajectories[i-1];

            io::Trajectory::const_iterator cit_b = std::lower_bound(trajectory.cbegin(), trajectory.cend(), timestamp, cmp);
            if (cit_b == trajectory.cend()) break;

            // Interpolate
            io::pose_t interpolated = sample(trajectory, timestamp);
            io::pose_t incremetal = relative_pose(last_poses[i-1], interpolated);

            p.push_back(incremetal);
            last_poses[i-1] = interpolated;
        }

        has_next = p.size() == incremental_poses.size();
        if (has_next) {
            // Save incremental motions
            for (int i = 1; i < argc; ++i)
                incremental_poses[i-1].push_back(p[i-1]);
            k = next_k;
        }
    } while (has_next);

    // Save incremental motions to file
    boost::filesystem::path output_path = FLAGS_output_dir.empty() ? boost::filesystem::current_path() : FLAGS_output_dir;
    for (int i = 1; i < argc; ++i) {
        boost::filesystem::path file_path = output_path;
        file_path /= (std::to_string(i) + ".txt");

        CHECK(!boost::filesystem::exists(file_path)) << "File already exists:" << std::endl
                                                     << file_path.string();
        CHECK(io::write_file<io::pose_t>(incremental_poses[i-1], file_path.string())) << "Unable to save output file:" << std::endl
                                                                                      << file_path.string();
    }

    return 0;
}
