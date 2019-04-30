#ifndef IO_HPP_
#define IO_HPP_

// STL
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace io {

// Types

using record_t = std::pair<double, std::string>;

struct pose_t {
    double tx, ty, theta;

    pose_t()
        : tx(0.0), ty(0.0), theta(0.0)
    { }

    pose_t(double tx, double ty, double theta)
        : tx(tx), ty(ty), theta(theta)
    { }
};

template<typename T>
struct trajectory_t {
    T id;
    pose_t pose;
};

using Records = std::vector<record_t>;
using Trajectory = std::vector<trajectory_t<double>>;

// ------------------------------------------------------------

// Comparators

template<typename T>
inline bool operator<(const trajectory_t<T> &lhs, const trajectory_t<T> &rhs) {
    return (lhs.id < rhs.id);
}

// ------------------------------------------------------------

// IO helpers

inline std::string to_string(int n, int w) {

    std::stringstream ss;
    ss << std::setw(w) << std::setfill('0') << n;
    return ss.str();
}

inline std::string to_string(double n, int precision) {

    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << n;
    return ss.str();
}

inline std::istream& operator>>(std::istream &lhs, record_t &rhs) {

    lhs >> rhs.first >> rhs.second;
    return lhs;
}

inline std::istream& operator>>(std::istream &lhs, pose_t &rhs) {

    lhs >> rhs.tx >> rhs.ty >> rhs.theta;
    return lhs;
}

template<typename T>
inline std::istream& operator>>(std::istream &lhs, trajectory_t<T> &rhs) {

    lhs >> rhs.id >> rhs.pose;
    return lhs;
}

inline std::ostream& operator<<(std::ostream &lhs, const record_t &rhs) {

    lhs << to_string(rhs.first, 9) << " " << rhs.second;
    return lhs;
}

inline std::ostream& operator<<(std::ostream &lhs, const pose_t &rhs) {

    lhs << rhs.tx << " " << rhs.ty << " " << rhs.theta;
    return lhs;
}

inline std::ostream& operator<<(std::ostream &lhs, const trajectory_t<double> &rhs) {

    lhs << to_string(rhs.id, 9) << " " << rhs.pose;
    return lhs;
}

inline std::ostream& operator<<(std::ostream &lhs, const trajectory_t<int> &rhs) {

    lhs << rhs.id << " " << rhs.pose;
    return lhs;
}

template<typename T>
inline std::vector<T> read_file(const std::string &path) {

    std::vector<T> records;

    std::ifstream input(path);
    if (!input.is_open()) return records;

    for (std::string line; std::getline(input, line);) {
        if (line.empty() || line.front() == '#') continue;

        std::istringstream iss(line);
        T record;
        if (iss >> record) records.push_back(std::move(record));
    }

    std::sort(records.begin(), records.end());
    return records;
}

template<typename T>
inline void write_stream(const std::vector<T> &records, std::ostream &stream) {
    for (const T &record : records)
        stream << record << std::endl;
}

template<typename T>
inline bool write_file(const std::vector<T> &records, const std::string &path) {

    std::ofstream output(path);
    if (!output.is_open()) return false;

    write_stream(records, output);

    output.close();
    return (!output.fail() && !output.bad());
}

// ------------------------------------------------------------

} // namespace io

#endif // IO_HPP_
