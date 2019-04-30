#ifndef IO_3D_HPP_
#define IO_3D_HPP_

#include "util/io.hpp"

#include <Eigen/Geometry>

namespace io {

struct pose3d_t {
    double tx, ty, tz, qw, qx, qy, qz;

    pose3d_t()
        : tx(0.0), ty(0.0), tz(0.0), qw(1.0), qx(0.0), qy(0.0), qz(0.0)
    { }

    pose3d_t(double tx, double ty, double tz, double qw, double qx, double qy, double qz)
        : tx(tx), ty(ty), tz(tz), qw(qw), qx(qx), qy(qy), qz(qz)
    {
        q_normalize();
    }

    // From Eigen Transform (constructor)
    template<typename Scalar, int Type>
    pose3d_t(const Eigen::Transform<Scalar, 3, Type> &T) {
        Eigen::Quaternion<Scalar> q(T.rotation());
        q.normalize();

        tx = T.translation()(0);
        ty = T.translation()(1);
        tz = T.translation()(2);
        qw = q.w();
        qx = q.x();
        qy = q.y();
        qz = q.z();
    }

    // To Eigen Transform (implicit conversion)
    template<typename Scalar, int Type>
    operator Eigen::Transform<Scalar, 3, Type>() const {
        Eigen::Quaternion<Scalar> q(qw, qx, qy, qz);
        q.normalize();

        Eigen::Transform<Scalar, 3, Type> T(q);

        T.translation()(0) = tx;
        T.translation()(1) = ty;
        T.translation()(2) = tz;

        return T;
    }

    // Inverse
    inline pose3d_t inverse() {
        return Eigen::Isometry3d(*this).inverse();
    }

private:
    inline void q_normalize() {
        double norm = std::sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
        qw /= norm;
        qx /= norm;
        qy /= norm;
        qz /= norm;
    }
};

template<typename T>
struct trajectory3d_t {
    T id;
    pose3d_t pose;

    trajectory3d_t()
        : id(0), pose()
    { }

    trajectory3d_t(T id, const pose3d_t& pose)
        : id(id), pose(pose)
    { }
};

using Trajectory3d = std::vector<trajectory3d_t<timestamp_t>>;

template<typename T>
inline bool operator<(const trajectory3d_t<T>& lhs, const trajectory3d_t<T>& rhs) {
    return (lhs.id < rhs.id);
}

inline std::ostream& operator<<(std::ostream& lhs, const pose3d_t& rhs) {

    lhs << to_string(rhs.tx, 9) << " " << to_string(rhs.ty, 9) << " " << to_string(rhs.tz, 9) << " "
        << to_string(rhs.qx, 9) << " " << to_string(rhs.qy, 9) << " " << to_string(rhs.qz, 9) << " " << to_string(rhs.qw, 9);
    return lhs;
}

inline std::ostream& operator<<(std::ostream &lhs, const trajectory3d_t<double> &rhs) {

    lhs << to_string(rhs.id, 9) << " " << rhs.pose;
    return lhs;
}

inline std::ostream& operator<<(std::ostream &lhs, const trajectory3d_t<int> &rhs) {

    lhs << rhs.id << " " << rhs.pose;
    return lhs;
}

} // namespace io

#endif // IO_3D_HPP_
