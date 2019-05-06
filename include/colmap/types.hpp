// Lightweight version of COLMAP basic types

#ifndef COLMAP_TYPES_HPP_
#define COLMAP_TYPES_HPP_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

// Eigen
#include <Eigen/Core>

#include "colmap/alignment.hpp"
#include "colmap/camera_models.h"
#include "colmap/math.hpp"
#include "colmap/pose.hpp"

namespace Eigen {

typedef Eigen::Matrix<double, 3, 4> Matrix3x4d;
typedef Eigen::Matrix<uint8_t, 3, 1> Vector3ub;
typedef Eigen::Matrix<uint8_t, 4, 1> Vector4ub;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

}  // namespace Eigen

namespace colmap {

// Define non-copyable or non-movable classes.
#define NON_COPYABLE(class_name)          \
  class_name(class_name const&) = delete; \
  void operator=(class_name const& obj) = delete;
#define NON_MOVABLE(class_name) class_name(class_name&&) = delete;

// Unique identifier for cameras.
using camera_t = std::uint32_t;

// Unique identifier for images.
using image_t = std::uint32_t;

// Each image pair gets a unique ID, see `Database::ImagePairToPairId`.
using image_pair_t = std::uint64_t;

// Index per image, i.e. determines maximum number of 2D points per image.
using point2D_t = std::uint32_t;

// Unique identifier per added 3D point. Since we add many 3D points,
// delete them, and possibly re-add them again, the maximum number of allowed
// unique indices should be large.
using point3D_t = std::uint64_t;

// Values for invalid identifiers or indices.
const camera_t kInvalidCameraId = std::numeric_limits<camera_t>::max();
const image_t kInvalidImageId = std::numeric_limits<image_t>::max();

const point2D_t kInvalidPoint2DIdx = std::numeric_limits<point2D_t>::max();
const point3D_t kInvalidPoint3DId = std::numeric_limits<point3D_t>::max();

// 2D point class corresponds to a feature in an image. It may or may not have a
// corresponding 3D point if it is part of a triangulated track.
class Point2D {
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Point2D()
        : xy_(Eigen::Vector2d::Zero()), point3D_id_(kInvalidPoint3DId) { }

    // The coordinate in image space in pixels.
    inline const Eigen::Vector2d& XY() const { return xy_; }
    inline Eigen::Vector2d& XY() { return xy_; }
    inline double X() const { return xy_.x(); }
    inline double Y() const { return xy_.y(); }
    inline void SetXY(const Eigen::Vector2d& xy) { xy_ = xy; }

    // The identifier of the observed 3D point. If the image point does not
    // observe a 3D point, the identifier is `kInvalidPoint3Did`.
    inline point3D_t Point3DId() const { return point3D_id_; }
    inline bool HasPoint3D() const { return point3D_id_ != kInvalidPoint3DId; }
    inline void SetPoint3DId(const point3D_t point3D_id) { point3D_id_ = point3D_id; }

private:
    // The image coordinates in pixels, starting at upper left corner with 0.
    Eigen::Vector2d xy_;

    // The identifier of the 3D point. If the 2D point is not part of a 3D point
    // track the identifier is `kInvalidPoint3DId` and `HasPoint3D() = false`.
    point3D_t point3D_id_;
};

// Track class stores all observations of a 3D point.
struct TrackElement {
    TrackElement()
        : image_id(kInvalidImageId), point2D_idx(kInvalidPoint2DIdx) { }

    TrackElement(const image_t image_id, const point2D_t point2D_idx)
        : image_id(image_id), point2D_idx(point2D_idx) { }

    // The image in which the track element is observed.
    image_t image_id;

    // The point in the image that the track element is observed.
    point2D_t point2D_idx;
};

class Track {
public:
    Track() { }

    // The number of track elements.
    inline std::size_t Length() const { return elements_.size(); }

    // Access all elements.
    inline const std::vector<TrackElement>& Elements() const { return elements_; }
    inline void SetElements(const std::vector<TrackElement>& elements) { elements_ = elements; }

    // Access specific elements.
    inline const TrackElement& Element(const std::size_t idx) const { return elements_.at(idx); }
    inline TrackElement& Element(const std::size_t idx) { return elements_.at(idx); }
    inline void SetElement(const std::size_t idx, const TrackElement& element) { elements_.at(idx) = element; }

    // Append new elements.
    inline void AddElement(const TrackElement& element) { elements_.push_back(element); }
    inline void AddElement(const image_t image_id, const point2D_t point2D_idx) { elements_.emplace_back(image_id, point2D_idx); }
    inline void AddElements(const std::vector<TrackElement>& elements) { elements_.insert(elements_.end(), elements.begin(), elements.end()); }

    // Delete existing element.
    inline void DeleteElement(const std::size_t idx) { if (idx < elements_.size()) elements_.erase(elements_.begin() + idx); }
    inline void DeleteElement(const image_t image_id, const point2D_t point2D_idx) {
        elements_.erase(
            std::remove_if(elements_.begin(), elements_.end(),
                           [image_id, point2D_idx](const TrackElement& element) {
                             return element.image_id == image_id &&
                                    element.point2D_idx == point2D_idx;
                           }),
            elements_.end());
    }

    // Requests that the track capacity be at least enough to contain the
    // specified number of elements.
    inline void Reserve(const std::size_t num_elements) { elements_.reserve(num_elements); }

    // Shrink the capacity of track vector to fit its size to save memory.
    inline void Compress() { elements_.shrink_to_fit(); }

private:
    std::vector<TrackElement> elements_;
};

// 3D point class that holds information about triangulated 2D points.
class Point3D {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Point3D()
        : xyz_(0.0, 0.0, 0.0) { }

    // The point coordinate in world space.
    inline const Eigen::Vector3d& XYZ() const { return xyz_; }
    inline Eigen::Vector3d& XYZ() { return xyz_; }
    inline double XYZ(const std::size_t idx) const { return xyz_(idx); }
    inline double& XYZ(const std::size_t idx) { return xyz_(idx); }
    inline double X() const { return xyz_.x(); }
    inline double Y() const { return xyz_.y(); }
    inline double Z() const { return xyz_.z(); }
    inline void SetXYZ(const Eigen::Vector3d& xyz) { xyz_ = xyz; }

    // The RGB color of the point.
    inline const Eigen::Vector3ub& Color() const { return color_; }
    inline Eigen::Vector3ub& Color() { return color_; }
    inline uint8_t Color(const std::size_t idx) const { return color_(idx); }
    inline uint8_t& Color(const std::size_t idx) { return color_(idx); }
    inline void SetColor(const Eigen::Vector3ub& color) { color_ = color; }

    // The mean reprojection error in image space.
    inline double Error() const { return error_; }
    inline bool HasError() const { return error_ != -1.0; }
    inline void SetError(const double error) { error_ = error; }

    inline const class Track& Track() const { return track_; }
    inline class Track& Track() { return track_; }
    inline void SetTrack(const class Track& track) { track_ = track; }

private:
    // The 3D position of the point.
    Eigen::Vector3d xyz_;

    // The color of the point in the range [0, 255].
    Eigen::Vector3ub color_;

    // The mean reprojection error in pixels.
    double error_;

    // The track of the point as a list of image observations.
    class Track track_;
};

// Class that holds information about an image. An image is the product of one
// camera shot at a certain location (parameterized as the pose). An image may
// share a camera with multiple other images, if its intrinsics are the same.
class Image {
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Image()
        : image_id_(kInvalidImageId),
          name_(""),
          camera_id_(kInvalidCameraId),
          registered_(false),
          num_points3D_(0),
          qvec_(ComposeIdentityQuaternion()),
          tvec_(Eigen::Vector3d::Zero()) { }

    // Access the unique identifier of the image.
    inline image_t ImageId() const { return image_id_; }
    inline void SetImageId(const image_t image_id) { image_id_ = image_id; }

    // Access the name of the image.
    inline const std::string& Name() const { return name_; }
    inline std::string& Name() { return name_; }
    inline void SetName(const std::string& name) { name_ = name; }

    // Access the unique identifier of the camera. Note that multiple images
    // might share the same camera.
    inline camera_t CameraId() const { return camera_id_; }
    inline void SetCameraId(const camera_t camera_id) { camera_id_ = camera_id; }
    // Check whether identifier of camera has been set.
    inline bool HasCamera() const { return camera_id_ != kInvalidCameraId; }

    // Check if image is registered.
    inline bool IsRegistered() const { return registered_; }
    inline void SetRegistered(const bool registered) { registered_ = registered; }

    // Get the number of image points.
    inline point2D_t NumPoints2D() const { return static_cast<point2D_t>(points2D_.size()); }

    // Get the number of triangulations, i.e. the number of points that
    // are part of a 3D point track.
    inline point2D_t NumPoints3D() const { return num_points3D_; }

    // Access quaternion vector as (qw, qx, qy, qz) specifying the rotation of the
    // pose which is defined as the transformation from world to image space.
    inline const Eigen::Vector4d& Qvec() const { return qvec_; }
    inline Eigen::Vector4d& Qvec() { return qvec_; }
    inline double Qvec(const std::size_t idx) const { return qvec_(idx); }
    inline double& Qvec(const std::size_t idx) { return qvec_(idx); }
    inline void SetQvec(const Eigen::Vector4d& qvec) { qvec_ = qvec; }

    // Quaternion prior, e.g. given by EXIF gyroscope tag.
    inline const Eigen::Vector4d& QvecPrior() const { return qvec_prior_; }
    inline Eigen::Vector4d& QvecPrior() { return qvec_prior_; }
    inline double QvecPrior(const size_t idx) const { return qvec_prior_(idx); }
    inline double& QvecPrior(const size_t idx) { return qvec_prior_(idx); }
    inline bool HasQvecPrior() const { return !IsNaN(qvec_prior_.sum()); }
    inline void SetQvecPrior(const Eigen::Vector4d& qvec) { qvec_prior_ = qvec; }

    // Access quaternion vector as (tx, ty, tz) specifying the translation of the
    // pose which is defined as the transformation from world to image space.
    inline const Eigen::Vector3d& Tvec() const { return tvec_; }
    inline Eigen::Vector3d& Tvec() { return tvec_; }
    inline double Tvec(const std::size_t idx) const { return tvec_(idx); }
    inline double& Tvec(const std::size_t idx) { return tvec_(idx); }
    inline void SetTvec(const Eigen::Vector3d& tvec) { tvec_ = tvec; }

    // Translation prior, e.g. given by EXIF GPS tag.
    inline const Eigen::Vector3d& TvecPrior() const { return tvec_prior_; }
    inline Eigen::Vector3d& TvecPrior() { return tvec_prior_; }
    inline double TvecPrior(const size_t idx) const { return tvec_prior_(idx); }
    inline double& TvecPrior(const size_t idx) { return tvec_prior_(idx); }
    inline bool HasTvecPrior() const { return !IsNaN(tvec_prior_.sum()); }
    inline void SetTvecPrior(const Eigen::Vector3d& tvec) { tvec_prior_ = tvec; }

    // Access the coordinates of image points.
    inline const class Point2D& Point2D(const point2D_t point2D_idx) const { return points2D_.at(point2D_idx); }
    inline class Point2D& Point2D(const point2D_t point2D_idx) { return points2D_.at(point2D_idx); }
    inline const std::vector<class Point2D>& Points2D() const { return points2D_; }
    inline void SetPoints2D(const std::vector<Eigen::Vector2d>& points) {
        points2D_.resize(points.size());
        for (point2D_t point2D_idx = 0; point2D_idx < points.size(); ++point2D_idx)
            points2D_[point2D_idx].SetXY(points[point2D_idx]);
    }
    inline void SetPoints2D(const std::vector<class Point2D>& points) { points2D_ = points; }

    // Set the point as triangulated, i.e. it is part of a 3D point track.
    inline void SetPoint3DForPoint2D(const point2D_t point2D_idx, const point3D_t point3D_id) {
        if (point3D_id == kInvalidPoint3DId) throw std::domain_error("Invalid 3D point id");
        class Point2D& point2D = points2D_.at(point2D_idx);
        if (!point2D.HasPoint3D()) {
          num_points3D_ += 1;
        }
        point2D.SetPoint3DId(point3D_id);
    }

    // Normalize the quaternion vector.
    inline void NormalizeQvec() { qvec_ = NormalizeQuaternion(qvec_); }

    // ProjectionCenterFromParameters
    inline Eigen::Vector3d ProjectionCenter() const {
        return ProjectionCenterFromParameters(qvec_, tvec_);
    }

private:
    // Identifier of the image, if not specified `kInvalidImageId`.
    image_t image_id_;

    // The name of the image, i.e. the relative path.
    std::string name_;

    // The identifier of the associated camera. Note that multiple images might
    // share the same camera. If not specified `kInvalidCameraId`.
    camera_t camera_id_;

    // Whether the image is successfully registered in the reconstruction.
    bool registered_;

    // The number of 3D points the image observes, i.e. the sum of its `points2D`
    // where `point3D_id != kInvalidPoint3DId`.
    point2D_t num_points3D_;

    // The pose of the image, defined as the transformation from world to image.
    Eigen::Vector4d qvec_;
    Eigen::Vector3d tvec_;

    // The pose prior of the image, e.g. extracted from EXIF tags.
    Eigen::Vector4d qvec_prior_;
    Eigen::Vector3d tvec_prior_;

    // All image points, including points that are not part of a 3D point track.
    std::vector<class Point2D> points2D_;
};

class Camera {
public:
    Camera()
        : camera_id_(kInvalidCameraId),
          model_id_(kInvalidCameraModelId),
          width_(0),
          height_(0) { }

    // Access the unique identifier of the camera.
    inline camera_t CameraId() const { return camera_id_; }
    inline void SetCameraId(const camera_t camera_id) { camera_id_ = camera_id; }

    // Access the camera model.
    inline int ModelId() const { return model_id_; }
    inline std::string ModelName() const { return CameraModelIdToName(model_id_); }
    inline void SetModelId(const int model_id) {
        if (!ExistsCameraModelWithId(model_id)) CAMERA_MODEL_DOES_NOT_EXIST_EXCEPTION;
        model_id_ = model_id;
        params_.resize(CameraModelNumParams(model_id_), 0);
    }
    inline void SetModelIdFromName(const std::string& model_name) {
        if (!ExistsCameraModelWithName(model_name)) CAMERA_MODEL_DOES_NOT_EXIST_EXCEPTION;
        model_id_ = CameraModelNameToId(model_name);
        params_.resize(CameraModelNumParams(model_id_), 0);
    }

    // Access dimensions of the camera sensor.
    inline std::size_t Width() const { return width_; }
    inline std::size_t Height() const { return height_; }
    inline void SetWidth(const std::size_t width) { width_ = width; }
    inline void SetHeight(const std::size_t height) { height_ = height; }

    // Get the indices of the parameter groups in the parameter vector.
    inline const std::vector<std::size_t>& FocalLengthIdxs() const { return CameraModelFocalLengthIdxs(model_id_); }
    inline const std::vector<std::size_t>& PrincipalPointIdxs() const { return CameraModelPrincipalPointIdxs(model_id_); }
    inline const std::vector<std::size_t>& ExtraParamsIdxs() const { return CameraModelExtraParamsIdxs(model_id_); }

    // Get human-readable information about the parameter vector ordering.
    inline std::string ParamsInfo() const { return CameraModelParamsInfo(model_id_); }

    // Access the raw parameter vector.
    inline std::size_t NumParams() const { return params_.size(); }
    inline const std::vector<double>& Params() const { return params_; }
    inline std::vector<double>& Params() { return params_; }
    inline double Params(const std::size_t idx) const { return params_[idx]; }
    inline double& Params(const std::size_t idx) { return params_[idx]; }
    inline const double* ParamsData() const { return params_.data(); }
    inline double* ParamsData() { return params_.data(); }
    inline void SetParams(const std::vector<double>& params) { params_ = params; }

    // Check whether parameters are valid, i.e. the parameter vector has
    // the correct dimensions that match the specified camera model.
    inline bool VerifyParams() const { return CameraModelVerifyParams(model_id_, params_); }

    // Project point in image plane to world / infinity.
    inline Eigen::Vector2d ImageToWorld(const Eigen::Vector2d& image_point) const {
        Eigen::Vector2d world_point;
        CameraModelImageToWorld(model_id_, params_, image_point(0), image_point(1),
                                &world_point(0), &world_point(1));
        return world_point;
    }

    // Project point from world / infinity to image plane.
    inline Eigen::Vector2d WorldToImage(const Eigen::Vector2d& world_point) const {
        Eigen::Vector2d image_point;
        CameraModelWorldToImage(model_id_, params_, world_point(0), world_point(1),
                                &image_point(0), &image_point(1));
        return image_point;
    }

private:
    // The unique identifier of the camera. If the identifier is not specified
    // it is set to `kInvalidCameraId`.
    camera_t camera_id_;

    // The identifier of the camera model. If the camera model is not specified
    // the identifier is `kInvalidCameraModelId`.
    int model_id_;

    // The dimensions of the image, 0 if not initialized.
    std::size_t width_;
    std::size_t height_;

    // The focal length, principal point, and extra parameters. If the camera
    // model is not specified, this vector is empty.
    std::vector<double> params_;
};

struct FeatureKeypoint {
    FeatureKeypoint()
        : FeatureKeypoint(0, 0) { }
    FeatureKeypoint(const float x, const float y)
        : FeatureKeypoint(x, y, 1, 0, 0, 1) { }
    FeatureKeypoint(const float x_, const float y_, const float scale,
                    const float orientation)
        : x(x_), y(y_) {
        if (scale < 0.0) throw std::domain_error("Invalid scale");
        const float scale_cos_orientation = scale * std::cos(orientation);
        const float scale_sin_orientation = scale * std::sin(orientation);
        a11 = scale_cos_orientation;
        a12 = -scale_sin_orientation;
        a21 = scale_sin_orientation;
        a22 = scale_cos_orientation;
    }
    FeatureKeypoint(const float x_, const float y_, const float a11_,
                    const float a12_, const float a21_, const float a22_)
        : x(x_), y(y_), a11(a11_), a12(a12_), a21(a21_), a22(a22_) { }

//  static FeatureKeypoint FromParameters(const float x, const float y,
//                                        const float scale_x,
//                                        const float scale_y,
//                                        const float orientation,
//                                        const float shear);

//  // Rescale the feature location and shape size by the given scale factor.
//  void Rescale(const float scale);
//  void Rescale(const float scale_x, const float scale_y);

//  // Compute similarity shape parameters from affine shape.
//  float ComputeScale() const;
//  float ComputeScaleX() const;
//  float ComputeScaleY() const;
//  float ComputeOrientation() const;
//  float ComputeShear() const;

    // Location of the feature, with the origin at the upper left image corner,
    // i.e. the upper left pixel has the coordinate (0.5, 0.5).
    float x;
    float y;

    // Affine shape of the feature.
    float a11;
    float a12;
    float a21;
    float a22;
};

typedef Eigen::Matrix<uint8_t, 1, Eigen::Dynamic, Eigen::RowMajor> FeatureDescriptor;

struct FeatureMatch {
    // Feature index in first image.
    point2D_t point2D_idx1 = kInvalidPoint2DIdx;
    // Feature index in second image.
    point2D_t point2D_idx2 = kInvalidPoint2DIdx;
};

typedef std::vector<FeatureKeypoint> FeatureKeypoints;
typedef Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> FeatureDescriptors;
typedef std::vector<FeatureMatch> FeatureMatches;

}  // namespace colmap

#endif  // COLMAP_TYPES_HPP_
