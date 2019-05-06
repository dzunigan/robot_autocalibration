// Minimal version of COLMAP reconstruction class

#ifndef COLMAP_SRC_BASE_RECONSTRUCTION_H_
#define COLMAP_SRC_BASE_RECONSTRUCTION_H_

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

// Boost
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include "colmap/alignment.hpp"
#include "colmap/endian.h"
#include "colmap/misc.h"
#include "colmap/types.hpp"

namespace colmap {

class Reconstruction {
public:

    Reconstruction()
        : num_added_points3D_(0) { }

    // Get const objects.
    inline const class Camera& Camera(const camera_t camera_id) const { return cameras_.at(camera_id); }
    inline const class Image& Image(const image_t image_id) const { return images_.at(image_id); }
    inline const class Point3D& Point3D(const point3D_t point3D_id) const { return points3D_.at(point3D_id); }

    // Get mutable objects.
    inline class Camera& Camera(const camera_t camera_id) { return cameras_.at(camera_id); }
    inline class Image& Image(const image_t image_id) { return images_.at(image_id); }
    inline class Point3D& Point3D(const point3D_t point3D_id) { return points3D_.at(point3D_id); }

    // Get reference to all objects.
    inline const EIGEN_STL_UMAP(camera_t, class Camera) & Cameras() const { return cameras_; }
    inline const EIGEN_STL_UMAP(image_t, class Image) & Images() const { return images_; }
    inline const std::vector<image_t>& RegImageIds() const { return reg_image_ids_; }
    inline const EIGEN_STL_UMAP(point3D_t, class Point3D) & Points3D() const { return points3D_; }

    // Check whether specific object exists.
    inline bool ExistsCamera(const camera_t camera_id) const { return cameras_.find(camera_id) != cameras_.end(); }
    inline bool ExistsImage(const image_t image_id) const { return images_.find(image_id) != images_.end(); }
    inline bool ExistsPoint3D(const point3D_t point3D_id) const { return points3D_.find(point3D_id) != points3D_.end(); }

    // Add new camera. There is only one camera per image, while multiple images
    // might be taken by the same camera.
    inline void AddCamera(const class Camera& camera) {
        if (ExistsCamera(camera.CameraId())) throw std::domain_error("Camera already added");
        if (!camera.VerifyParams()) throw std::runtime_error("Invalid camera parameters");
        cameras_.emplace(camera.CameraId(), camera);
    }

    // Add new image.
    inline void AddImage(const class Image& image) {
        if (ExistsImage(image.ImageId())) throw std::domain_error("Image already added");
        images_[image.ImageId()] = image;
    }

    // Add new 3D object, and return its unique ID.
    inline point3D_t AddPoint3D(const Eigen::Vector3d& xyz, const Track& track) {
        const point3D_t point3D_id = ++num_added_points3D_;
        if (ExistsPoint3D(point3D_id)) throw std::runtime_error("Invalid 3D point id");

        class Point3D& point3D = points3D_[point3D_id];

        point3D.SetXYZ(xyz);
        point3D.SetTrack(track);

        for (const auto& track_el : track.Elements()) {
          class Image& image = Image(track_el.image_id);
          if (image.Point2D(track_el.point2D_idx).HasPoint3D()) throw std::domain_error("Keypoint has 3D information");
          image.SetPoint3DForPoint2D(track_el.point2D_idx, point3D_id);
          if (image.NumPoints3D() > image.NumPoints2D()) throw std::runtime_error("Invalid number of 3d points");
        }

//        const bool kIsContinuedPoint3D = false;

//        for (const auto& track_el : track.Elements()) {
//          SetObservationAsTriangulated(track_el.image_id, track_el.point2D_idx,
//                                       kIsContinuedPoint3D);
//        }
    }

    // Add observation to existing 3D point.
    inline void AddObservation(const point3D_t point3D_id, const TrackElement& track_el) {
        class Image& image = Image(track_el.image_id);
        if (image.Point2D(track_el.point2D_idx).HasPoint3D()) throw std::domain_error("Keypoint has 3D information");

        image.SetPoint3DForPoint2D(track_el.point2D_idx, point3D_id);
        if (image.NumPoints3D() > image.NumPoints2D()) throw std::runtime_error("Invalid number of 3d points");

        class Point3D& point3D = Point3D(point3D_id);
        point3D.Track().AddElement(track_el);

//        const bool kIsContinuedPoint3D = true;
//        SetObservationAsTriangulated(track_el.image_id, track_el.point2D_idx,
//                                     kIsContinuedPoint3D);
    }

    // Normalize scene by scaling and translation to avoid degenerate
    // visualization after bundle adjustment and to improve numerical
    // stability of algorithms.
    //
    // Translates scene such that the mean of the camera centers or point
    // locations are at the origin of the coordinate system.
    //
    // Scales scene such that the minimum and maximum camera centers are at the
    // given `extent`, whereas `p0` and `p1` determine the minimum and
    // maximum percentiles of the camera centers considered.
    inline void Normalize(const double extent = 10.0, const double p0 = 0.1,
                          const double p1 = 0.9, const bool use_images = true) {
        if (extent <= 0 ||
                p0 < 0 || p0 > 1 ||
                p1 < 0 || p1 > 1 ||
                p0 >= p1)
            throw std::domain_error("Invalid parameter");

        if ((use_images && reg_image_ids_.size() < 2) ||
            (!use_images && points3D_.size() < 2)) {
          return;
        }

        EIGEN_STL_UMAP(class Image*, Eigen::Vector3d) proj_centers;

        for (size_t i = 0; i < reg_image_ids_.size(); ++i) {
          class Image& image = Image(reg_image_ids_[i]);
          const Eigen::Vector3d proj_center = image.ProjectionCenter();
          proj_centers[&image] = proj_center;
        }

        // Coordinates of image centers or point locations.
        std::vector<float> coords_x;
        std::vector<float> coords_y;
        std::vector<float> coords_z;
        if (use_images) {
          coords_x.reserve(proj_centers.size());
          coords_y.reserve(proj_centers.size());
          coords_z.reserve(proj_centers.size());
          for (const auto& proj_center : proj_centers) {
            coords_x.push_back(static_cast<float>(proj_center.second(0)));
            coords_y.push_back(static_cast<float>(proj_center.second(1)));
            coords_z.push_back(static_cast<float>(proj_center.second(2)));
          }
        } else {
          coords_x.reserve(points3D_.size());
          coords_y.reserve(points3D_.size());
          coords_z.reserve(points3D_.size());
          for (const auto& point3D : points3D_) {
            coords_x.push_back(static_cast<float>(point3D.second.X()));
            coords_y.push_back(static_cast<float>(point3D.second.Y()));
            coords_z.push_back(static_cast<float>(point3D.second.Z()));
          }
        }

        // Determine robust bounding box and mean.

        std::sort(coords_x.begin(), coords_x.end());
        std::sort(coords_y.begin(), coords_y.end());
        std::sort(coords_z.begin(), coords_z.end());

        const size_t P0 = static_cast<size_t>(
            (coords_x.size() > 3) ? p0 * (coords_x.size() - 1) : 0);
        const size_t P1 = static_cast<size_t>(
            (coords_x.size() > 3) ? p1 * (coords_x.size() - 1) : coords_x.size() - 1);

        const Eigen::Vector3d bbox_min(coords_x[P0], coords_y[P0], coords_z[P0]);
        const Eigen::Vector3d bbox_max(coords_x[P1], coords_y[P1], coords_z[P1]);

        Eigen::Vector3d mean_coord(0, 0, 0);
        for (size_t i = P0; i <= P1; ++i) {
          mean_coord(0) += coords_x[i];
          mean_coord(1) += coords_y[i];
          mean_coord(2) += coords_z[i];
        }
        mean_coord /= P1 - P0 + 1;

        // Calculate scale and translation, such that
        // translation is applied before scaling.
        const double old_extent = (bbox_max - bbox_min).norm();
        double scale;
        if (old_extent < std::numeric_limits<double>::epsilon()) {
          scale = 1;
        } else {
          scale = extent / old_extent;
        }

        const Eigen::Vector3d translation = mean_coord;

        // Transform images.
        for (auto& elem : proj_centers) {
          elem.second -= translation;
          elem.second *= scale;
          const Eigen::Quaterniond quat(elem.first->Qvec(0), elem.first->Qvec(1),
                                        elem.first->Qvec(2), elem.first->Qvec(3));
          elem.first->SetTvec(quat * -elem.second);
        }

        // Transform points.
        for (auto& point3D : points3D_) {
          point3D.second.XYZ() -= translation;
          point3D.second.XYZ() *= scale;
        }
    }

    // Read data from text or binary file. Prefer binary data if it exists.
    inline void Read(const std::string& path) { ReadBinary(path); }
    inline void Write(const std::string& path) const { WriteBinary(path); }

    // Read data from binary file.
    inline void ReadBinary(const std::string& path) {
        ReadCamerasBinary(JoinPaths(path, "cameras.bin"));
        ReadImagesBinary(JoinPaths(path, "images.bin"));
        ReadPoints3DBinary(JoinPaths(path, "points3D.bin"));
    }

    // Write data from binary/text file.
//    void WriteText(const std::string& path) const;
    inline void WriteBinary(const std::string& path) const {
        WriteCamerasBinary(JoinPaths(path, "cameras.bin"));
        WriteImagesBinary(JoinPaths(path, "images.bin"));
        WritePoints3DBinary(JoinPaths(path, "points3D.bin"));
    }

private:

    inline void ReadCamerasBinary(const std::string& path) {
        std::ifstream file(path, std::ios::binary);
        if (!boost::filesystem::is_regular_file(path) || !file.is_open()) return;

        const size_t num_cameras = ReadBinaryLittleEndian<std::uint64_t>(&file);
        for (size_t i = 0; i < num_cameras; ++i) {
            class Camera camera;
            camera.SetCameraId(ReadBinaryLittleEndian<camera_t>(&file));
            camera.SetModelId(ReadBinaryLittleEndian<int>(&file));
            camera.SetWidth(ReadBinaryLittleEndian<uint64_t>(&file));
            camera.SetHeight(ReadBinaryLittleEndian<uint64_t>(&file));

            for (std::size_t i = 0; i < camera.NumParams(); ++i)
                camera.Params(i) = ReadBinaryLittleEndian<double>(&file);

            if (!camera.VerifyParams()) throw std::domain_error("Invalid number of camera parameters");
            cameras_.emplace(camera.CameraId(), camera);
        }
    }

    inline void ReadImagesBinary(const std::string& path) {
        std::ifstream file(path, std::ios::binary);
        if (!boost::filesystem::is_regular_file(path) || !file.is_open()) return;

        const std::size_t num_reg_images = ReadBinaryLittleEndian<std::uint64_t>(&file);
        for (std::size_t i = 0; i < num_reg_images; ++i) {
            class Image image;

            image.SetImageId(ReadBinaryLittleEndian<image_t>(&file));

            image.Qvec(0) = ReadBinaryLittleEndian<double>(&file);
            image.Qvec(1) = ReadBinaryLittleEndian<double>(&file);
            image.Qvec(2) = ReadBinaryLittleEndian<double>(&file);
            image.Qvec(3) = ReadBinaryLittleEndian<double>(&file);
            image.NormalizeQvec();

            image.Tvec(0) = ReadBinaryLittleEndian<double>(&file);
            image.Tvec(1) = ReadBinaryLittleEndian<double>(&file);
            image.Tvec(2) = ReadBinaryLittleEndian<double>(&file);

            image.SetCameraId(ReadBinaryLittleEndian<camera_t>(&file));

            for(char name_char;;) {
                file.read(&name_char, 1);
                if (name_char != '\0') image.Name() += name_char;
                else break;
            }

            const std::size_t num_points2D = ReadBinaryLittleEndian<std::uint64_t>(&file);

            std::vector<Eigen::Vector2d> points2D;
            points2D.reserve(num_points2D);
            std::vector<point3D_t> point3D_ids;
            point3D_ids.reserve(num_points2D);
            for (size_t j = 0; j < num_points2D; ++j) {
                const double x = ReadBinaryLittleEndian<double>(&file);
                const double y = ReadBinaryLittleEndian<double>(&file);
                points2D.emplace_back(x, y);
                point3D_ids.push_back(ReadBinaryLittleEndian<point3D_t>(&file));
            }

            image.SetPoints2D(points2D);

            for (point2D_t point2D_idx = 0; point2D_idx < image.NumPoints2D(); ++point2D_idx)
                if (point3D_ids[point2D_idx] != kInvalidPoint3DId)
                    image.SetPoint3DForPoint2D(point2D_idx, point3D_ids[point2D_idx]);

            image.SetRegistered(true);
            reg_image_ids_.push_back(image.ImageId());

            images_.emplace(image.ImageId(), image);
        }
    }

    inline void ReadPoints3DBinary(const std::string& path) {
        std::ifstream file(path, std::ios::binary);
        if (!boost::filesystem::is_regular_file(path) || !file.is_open()) return;

        const size_t num_points3D = ReadBinaryLittleEndian<uint64_t>(&file);
        for (size_t i = 0; i < num_points3D; ++i) {
            class Point3D point3D;

            const point3D_t point3D_id = ReadBinaryLittleEndian<point3D_t>(&file);
            num_added_points3D_ = std::max(num_added_points3D_, point3D_id);

            point3D.XYZ()(0) = ReadBinaryLittleEndian<double>(&file);
            point3D.XYZ()(1) = ReadBinaryLittleEndian<double>(&file);
            point3D.XYZ()(2) = ReadBinaryLittleEndian<double>(&file);
            point3D.Color(0) = ReadBinaryLittleEndian<uint8_t>(&file);
            point3D.Color(1) = ReadBinaryLittleEndian<uint8_t>(&file);
            point3D.Color(2) = ReadBinaryLittleEndian<uint8_t>(&file);
            point3D.SetError(ReadBinaryLittleEndian<double>(&file));

            const size_t track_length = ReadBinaryLittleEndian<uint64_t>(&file);
            for (size_t j = 0; j < track_length; ++j) {
                const image_t image_id = ReadBinaryLittleEndian<image_t>(&file);
                const point2D_t point2D_idx = ReadBinaryLittleEndian<point2D_t>(&file);
                point3D.Track().AddElement(image_id, point2D_idx);
            }
            point3D.Track().Compress();

            points3D_.emplace(point3D_id, point3D);
        }
    }

//    void WriteCamerasText(const std::string& path) const;
//    void WriteImagesText(const std::string& path) const;
//    void WritePoints3DText(const std::string& path) const;
    inline void WriteCamerasBinary(const std::string& path) const {
        std::ofstream file(path, std::ios::trunc | std::ios::binary);
        if (!file.is_open()) throw std::runtime_error("Unable to open output file: " + path);

        WriteBinaryLittleEndian<uint64_t>(&file, cameras_.size());

        for (const auto& camera : cameras_) {
            WriteBinaryLittleEndian<camera_t>(&file, camera.first);
            WriteBinaryLittleEndian<int>(&file, camera.second.ModelId());
            WriteBinaryLittleEndian<uint64_t>(&file, camera.second.Width());
            WriteBinaryLittleEndian<uint64_t>(&file, camera.second.Height());
            for (const double param : camera.second.Params())
                WriteBinaryLittleEndian<double>(&file, param);
        }
    }

    inline void WriteImagesBinary(const std::string& path) const {
        std::ofstream file(path, std::ios::trunc | std::ios::binary);
        if (!file.is_open()) throw std::runtime_error("Unable to open output file: " + path);

        WriteBinaryLittleEndian<uint64_t>(&file, reg_image_ids_.size());

        for (const auto& image : images_) {
            if (!image.second.IsRegistered()) continue;

            WriteBinaryLittleEndian<image_t>(&file, image.first);

            const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(image.second.Qvec());
            WriteBinaryLittleEndian<double>(&file, normalized_qvec(0));
            WriteBinaryLittleEndian<double>(&file, normalized_qvec(1));
            WriteBinaryLittleEndian<double>(&file, normalized_qvec(2));
            WriteBinaryLittleEndian<double>(&file, normalized_qvec(3));

            WriteBinaryLittleEndian<double>(&file, image.second.Tvec(0));
            WriteBinaryLittleEndian<double>(&file, image.second.Tvec(1));
            WriteBinaryLittleEndian<double>(&file, image.second.Tvec(2));

            WriteBinaryLittleEndian<camera_t>(&file, image.second.CameraId());

            const std::string name = image.second.Name() + '\0';
            file.write(name.c_str(), name.size());

            WriteBinaryLittleEndian<uint64_t>(&file, image.second.NumPoints2D());
            for (const Point2D& point2D : image.second.Points2D()) {
                WriteBinaryLittleEndian<double>(&file, point2D.X());
                WriteBinaryLittleEndian<double>(&file, point2D.Y());
                WriteBinaryLittleEndian<point3D_t>(&file, point2D.Point3DId());
            }
        }
    }

    inline void WritePoints3DBinary(const std::string& path) const {
        std::ofstream file(path, std::ios::trunc | std::ios::binary);
        if (!file.is_open()) throw std::runtime_error("Unable to open output file: " + path);

        WriteBinaryLittleEndian<uint64_t>(&file, points3D_.size());

        for (const auto& point3D : points3D_) {
            WriteBinaryLittleEndian<point3D_t>(&file, point3D.first);
            WriteBinaryLittleEndian<double>(&file, point3D.second.XYZ()(0));
            WriteBinaryLittleEndian<double>(&file, point3D.second.XYZ()(1));
            WriteBinaryLittleEndian<double>(&file, point3D.second.XYZ()(2));
            WriteBinaryLittleEndian<uint8_t>(&file, point3D.second.Color(0));
            WriteBinaryLittleEndian<uint8_t>(&file, point3D.second.Color(1));
            WriteBinaryLittleEndian<uint8_t>(&file, point3D.second.Color(2));
            WriteBinaryLittleEndian<double>(&file, point3D.second.Error());

            WriteBinaryLittleEndian<uint64_t>(&file, point3D.second.Track().Length());
            for (const auto& track_el : point3D.second.Track().Elements()) {
                WriteBinaryLittleEndian<image_t>(&file, track_el.image_id);
                WriteBinaryLittleEndian<point2D_t>(&file, track_el.point2D_idx);
            }
        }
    }

    EIGEN_STL_UMAP(camera_t, class Camera) cameras_;
    EIGEN_STL_UMAP(image_t, class Image) images_;
    EIGEN_STL_UMAP(point3D_t, class Point3D) points3D_;

    // { image_id, ... } where `images_.at(image_id).registered == true`.
    std::vector<image_t> reg_image_ids_;

    // Total number of added 3D points, used to generate unique identifiers.
    point3D_t num_added_points3D_;
};

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_RECONSTRUCTION_H_
