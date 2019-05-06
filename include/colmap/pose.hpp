// Lightweight version of COLMAP 3D pose representation

#ifndef COLMAP_SRC_BASE_POSE_H_
#define COLMAP_SRC_BASE_POSE_H_

// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace colmap {

// Compose the Quaternion vector corresponding to a  identity transformation.
inline Eigen::Vector4d ComposeIdentityQuaternion() {
    return Eigen::Vector4d(1, 0, 0, 0);
}

// Normalize Quaternion vector.
//
// @param qvec          Quaternion rotation coefficients (w, x, y, z).
//
// @return              Unit Quaternion rotation coefficients (w, x, y, z).
inline Eigen::Vector4d NormalizeQuaternion(const Eigen::Vector4d& qvec) {
    const double norm = qvec.norm();
    if (norm == 0) {
        // We do not just use (1, 0, 0, 0) because that is a constant and when used
        // for automatic differentiation that would lead to a zero derivative.
        return Eigen::Vector4d(1.0, qvec(1), qvec(2), qvec(3));
    } else {
        const double inv_norm = 1.0 / norm;
        return inv_norm * qvec;
    }
}

// Convert 3D rotation matrix to Quaternion representation.
//
// @param rot_mat        3x3 rotation matrix.
//
// @return               Unit Quaternion rotation coefficients (w, x, y, z).
inline Eigen::Vector4d RotationMatrixToQuaternion(const Eigen::Matrix3d& rot_mat) {
    const Eigen::Quaterniond quat(rot_mat);
    return Eigen::Vector4d(quat.w(), quat.x(), quat.y(), quat.z());
}

// Convert Quaternion representation to 3D rotation matrix.
//
// @param qvec           Unit Quaternion rotation coefficients (w, x, y, z).
//
// @return               3x3 rotation matrix.
inline Eigen::Matrix3d QuaternionToRotationMatrix(const Eigen::Vector4d& qvec) {
    const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
    const Eigen::Quaterniond quat(normalized_qvec(0), normalized_qvec(1),
                                  normalized_qvec(2), normalized_qvec(3));
    return quat.toRotationMatrix();
}

// Invert Quaternion vector to return Quaternion of inverse rotation.
//
// @param qvec          Quaternion rotation coefficients (w, x, y, z).
//
// @return              Inverse Quaternion rotation coefficients (w, x, y, z).
inline Eigen::Vector4d InvertQuaternion(const Eigen::Vector4d& qvec) {
    const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
    return Eigen::Vector4d(normalized_qvec(0), -normalized_qvec(1),
                           -normalized_qvec(2), -normalized_qvec(3));
}

// Concatenate Quaternion rotations such that the rotation of `qvec1` is applied
// before the rotation of `qvec2`.
//
// @param qvec1         Quaternion rotation coefficients (w, x, y, z).
// @param qvec2         Quaternion rotation coefficients (w, x, y, z).
//
// @return              Concatenated Quaternion coefficients (w, x, y, z).
inline Eigen::Vector4d ConcatenateQuaternions(const Eigen::Vector4d& qvec1,
                                              const Eigen::Vector4d& qvec2) {
    const Eigen::Vector4d normalized_qvec1 = NormalizeQuaternion(qvec1);
    const Eigen::Vector4d normalized_qvec2 = NormalizeQuaternion(qvec2);
    const Eigen::Quaterniond quat1(normalized_qvec1(0), normalized_qvec1(1),
                                   normalized_qvec1(2), normalized_qvec1(3));
    const Eigen::Quaterniond quat2(normalized_qvec2(0), normalized_qvec2(1),
                                   normalized_qvec2(2), normalized_qvec2(3));
    const Eigen::Quaterniond cat_quat = quat2 * quat1;
    return NormalizeQuaternion(
                Eigen::Vector4d(cat_quat.w(), cat_quat.x(), cat_quat.y(), cat_quat.z()));
}

// Transform point by quaternion rotation.
//
// @param qvec          Quaternion rotation coefficients (w, x, y, z).
// @param point         Point to rotate.
//
// @return              Rotated point.
inline Eigen::Vector3d QuaternionRotatePoint(const Eigen::Vector4d& qvec,
                                             const Eigen::Vector3d& point) {
    const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
    const Eigen::Quaterniond quat(normalized_qvec(0), normalized_qvec(1),
                                normalized_qvec(2), normalized_qvec(3));
    return quat * point;
}

// Compute the relative transformation from pose 1 to 2.
//
// @param qvec1, tvec1      First camera pose.
// @param qvec2, tvec2      Second camera pose.
// @param qvec12, tvec12    Relative pose.
inline void ComputeRelativePose(const Eigen::Vector4d& qvec1,
                                const Eigen::Vector3d& tvec1,
                                const Eigen::Vector4d& qvec2,
                                const Eigen::Vector3d& tvec2, Eigen::Vector4d* qvec12,
                                Eigen::Vector3d* tvec12) {
    const Eigen::Vector4d inv_qvec1 = InvertQuaternion(qvec1);
    *qvec12 = ConcatenateQuaternions(inv_qvec1, qvec2);
    *tvec12 = tvec2 - QuaternionRotatePoint(*qvec12, tvec1);
}

// Concatenate the transformations of the two poses.
//
// @param qvec1, tvec1      First camera pose.
// @param qvec2, tvec2      Second camera pose.
// @param qvec12, tvec12    Concatenated pose.
inline void ConcatenatePoses(const Eigen::Vector4d& qvec1,
                             const Eigen::Vector3d& tvec1,
                             const Eigen::Vector4d& qvec2,
                             const Eigen::Vector3d& tvec2, Eigen::Vector4d* qvec12,
                             Eigen::Vector3d* tvec12) {
    *qvec12 = ConcatenateQuaternions(qvec1, qvec2);
    *tvec12 = tvec2 + QuaternionRotatePoint(qvec2, tvec1);
}

// Invert transformation of the pose.
// @param qvec, tvec          Input camera pose.
// @param inv_qvec, inv_tvec  Inverse camera pose.
inline void InvertPose(const Eigen::Vector4d& qvec, const Eigen::Vector3d& tvec,
                       Eigen::Vector4d* inv_qvec, Eigen::Vector3d* inv_tvec) {
    *inv_qvec = InvertQuaternion(qvec);
    *inv_tvec = -QuaternionRotatePoint(*inv_qvec, tvec);
}

// Linearly interpolate camera pose.
//
// @param qvec1, tvec1      Camera pose at t0 = 0.
// @param qvec2, tvec2      Camera pose at t1 = 1.
// @param t                 Interpolation time.
// @param qveci, tveci      Camera pose at time t.
inline void InterpolatePose(const Eigen::Vector4d& qvec1, const Eigen::Vector3d& tvec1,
                            const Eigen::Vector4d& qvec2, const Eigen::Vector3d& tvec2,
                            const double t, Eigen::Vector4d* qveci,
                            Eigen::Vector3d* tveci) {
    const Eigen::Vector4d normalized_qvec1 = NormalizeQuaternion(qvec1);
    const Eigen::Vector4d normalized_qvec2 = NormalizeQuaternion(qvec2);
    const Eigen::Quaterniond quat1(normalized_qvec1(0), normalized_qvec1(1),
                                 normalized_qvec1(2), normalized_qvec1(3));
    const Eigen::Quaterniond quat2(normalized_qvec2(0), normalized_qvec2(1),
                                 normalized_qvec2(2), normalized_qvec2(3));
    const Eigen::Vector3d tvec12 = tvec2 - tvec1;

    const Eigen::Quaterniond quati = quat1.slerp(t, quat2);

    *qveci = Eigen::Vector4d(quati.w(), quati.x(), quati.y(), quati.z());
    *tveci = tvec1 + tvec12 * t;
}

// Extract camera projection center from projection parameters.
//
// @param qvec           Unit Quaternion rotation coefficients (w, x, y, z).
// @param tvec           3x1 translation vector.
//
// @return               3x1 camera projection center.
Eigen::Vector3d ProjectionCenterFromParameters(const Eigen::Vector4d& qvec,
                                               const Eigen::Vector3d& tvec) {
  // Inverse rotation as conjugate quaternion.
  const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
  const Eigen::Quaterniond quat(normalized_qvec(0), -normalized_qvec(1),
                                -normalized_qvec(2), -normalized_qvec(3));
  return quat * -tvec;
}

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_POSE_H_
