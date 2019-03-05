#pragma once

// STL
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

// Ceres
#include <ceres/ceres.h>

#include "util/alignment.h"
#include "util/logging.h"

#include "angle_local_parameterization.h"
#include "sim2.hpp"

// Types
using sensor_t = std::uint32_t;
using sensor_pait_t = std::uint64_t;

// Constants
constexpr sensor_t kMaxNumSensors = std::numeric_limits<sensor_t>::max();
constexpr sensor_t kInvalidSensorId = kMaxNumSensors;
constexpr sensor_pait_t kInvalidSensorPairId = std::numeric_limits<sensor_pait_t>::max();

// Basic calibration residual.
struct CalibrationResidual {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CalibrationResidual(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
        : a(a(0), a(1), a(2), 1.0), b(b(0), b(1), b(2), 1.0) { }

    template <typename T>
    bool operator()(const T* const x_ptr, const T* const y_ptr, const T* const theta_ptr, const T* const s_ptr,
                    T* residuals_ptr) const {
        // Calibration parameters map
        Eigen::Matrix<T, 4, 1> x(*x_ptr, *y_ptr, *theta_ptr, *s_ptr);

        // Compute the error vector
        //Eigen::Matrix<T, 4, 1> e = PoseComposition<T>(b.cast<T>(), PoseInverse<T>(x)) - PoseComposition<T>(PoseInverse<T>(x), a.cast<T>());
        Eigen::Matrix<T, 4, 1> e = a.cast<T>() - PoseComposition<T>(x, PoseComposition<T>(b.cast<T>(), PoseInverse<T>(x)));

        // Compute the residuals
        Eigen::Map<Eigen::Matrix<T, 2, 1> > residuals(residuals_ptr);
        residuals = e.template head<2>();

        return true;
    }

    static ceres::CostFunction* Create(const Eigen::Vector3d &a, const Eigen::Vector3d &b) {
        return (new ceres::AutoDiffCostFunction<CalibrationResidual, 2, 1, 1, 1, 1>(
                    new CalibrationResidual(a, b)));
    }

private:

    const Eigen::Vector4d a, b; // Observations from first and second frames, respectively
};

// Additional calibration residual.
struct AdditionalCalibrationResidual {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AdditionalCalibrationResidual(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
        : a(a(0), a(1), a(2), 1.0), b(b(0), b(1), b(2), 1.0) { }

    template <typename T>
    bool operator()(const T* const x1_ptr, const T* const y1_ptr, const T* const theta1_ptr, const T* const s1_ptr,
                    const T* const x2_ptr, const T* const y2_ptr, const T* const theta2_ptr, const T* const s2_ptr,
                    T* residuals_ptr) const {
        // Calibration parameters map
        Eigen::Matrix<T, 4, 1> x(*x1_ptr, *y1_ptr, *theta1_ptr, *s1_ptr);
        Eigen::Matrix<T, 4, 1> y(*x2_ptr, *y2_ptr, *theta2_ptr, *s2_ptr);

        // Compute the effective calibration parameters
        Eigen::Matrix<T, 4, 1> z = PoseComposition<T>(PoseInverse<T>(x), y);

        // Compute the error vector
        //Eigen::Matrix<T, 4, 1> e = PoseComposition<T>(b.cast<T>(), PoseInverse<T>(z)) - PoseComposition<T>(PoseInverse<T>(z), a.cast<T>());
        Eigen::Matrix<T, 4, 1> e = a.cast<T>() - PoseComposition<T>(z, PoseComposition<T>(b.cast<T>(), PoseInverse<T>(z)));

        // Compute the residuals
        Eigen::Map<Eigen::Matrix<T, 2, 1> > residuals(residuals_ptr);
        residuals = e.template head<2>();

        return true;
    }

    static ceres::CostFunction* Create(const Eigen::Vector3d &a, const Eigen::Vector3d &b) {
        return (new ceres::AutoDiffCostFunction<AdditionalCalibrationResidual, 2, 1, 1, 1, 1, 1, 1, 1, 1>(
                    new AdditionalCalibrationResidual(a, b)));
    }

private:

    const Eigen::Vector4d a, b; // Observations from first and second frames, respectively
};

// Container for sensor-related information
class Sensor {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Sensor()
        : Sensor(kInvalidSensorId) { }

    explicit Sensor(sensor_t sensor_id)
        : sensor_id_(sensor_id), calibration_(0, 0, 0, 1) { }

    // Access the unique identifier of the sensor
    inline sensor_t SensorId() const { return sensor_id_; }
    inline void SetSensorId(const sensor_t sensor_id) { sensor_id_ = sensor_id; }

    // Access calibration parameters as (x, y, theta, s) specifying the sensor
    // to common frame coordinates transformation
    inline const Eigen::Vector4d& Calibration() const { return calibration_; }
    inline Eigen::Vector4d& Calibration() { return calibration_; }
    inline double Calibration(const size_t idx) const { return calibration_(idx); }
    inline double& Calibration(const size_t idx) { return calibration_(idx); }
    inline void SetCalibration(const Eigen::Vector4d& calibration) { calibration_ = calibration; }

    // Access the incremental motion observed by the sensor as (dx, dy, dtheta)
    inline const std::vector<Eigen::Vector3d>& Observations() const { return observations_; }
    inline std::vector<Eigen::Vector3d>& Observations() { return observations_; }
    inline const Eigen::Vector3d& Obsevation(const size_t idx) const { return observations_.at(idx); }
    inline Eigen::Vector3d& Obsevation(const size_t idx) { return observations_.at(idx); }
    inline void SetObservations(const std::vector<Eigen::Vector3d>& observations) { observations_ = observations; }

private:
    sensor_t sensor_id_;
    Eigen::Vector4d calibration_;
    std::vector<Eigen::Vector3d> observations_;
};

// Index pair container
struct index_pair_t {

    size_t idx1, idx2;

    index_pair_t(const size_t idx1, const size_t idx2)
        : idx1(idx1), idx2(idx2) { }

    bool operator==(const index_pair_t& rhs) const {
        return (idx1 == rhs.idx1 && idx2 == rhs.idx2);
    }
};

// Specialization of std::hash with index_pair_t
namespace std {
    template <>
    struct hash<index_pair_t> {
        size_t operator()( const index_pair_t& key ) const {
            // https://stackoverflow.com/a/1646913
            size_t hash = 17;
            hash = hash * 31 + key.idx1;
            hash = hash * 31 + key.idx2;
            return hash;
        }
    };
}

// Container manager for corresponding incremental motions between sensors
class ObservationDatabase {
public:

    inline static bool SwapSensorPair(const sensor_t sensor_id1, const sensor_t sensor_id2) {
        return sensor_id1 > sensor_id2;
    }

    inline static sensor_pait_t SensorPairToPairId(const sensor_t sensor_id1, const sensor_t sensor_id2) {
        CHECK_GE(sensor_id1, 0);
        CHECK_GE(sensor_id2, 0);
        CHECK_LT(sensor_id1, kMaxNumSensors);
        CHECK_LT(sensor_id2, kMaxNumSensors);

        if (SwapSensorPair(sensor_id1, sensor_id2)) {
            return kMaxNumSensors * sensor_id2 + sensor_id1;
        } else {
            return kMaxNumSensors * sensor_id1 + sensor_id2;
        }
    }

    inline static void PairIdToSensorPair(const sensor_pait_t pair_id, sensor_t* sensor_id1, sensor_t* sensor_id2) {
        *sensor_id2 = static_cast<sensor_t>(pair_id % kMaxNumSensors);
        *sensor_id1 = static_cast<sensor_t>((pair_id - *sensor_id2) / kMaxNumSensors);

        CHECK_GE(*sensor_id1, 0);
        CHECK_GE(*sensor_id2, 0);
        CHECK_LT(*sensor_id1, kMaxNumSensors);
        CHECK_LT(*sensor_id2, kMaxNumSensors);
    }

    ObservationDatabase() { }

    void AddIndexPair(const sensor_t sensor_id1, const sensor_t sensor_id2, const size_t idx1, const size_t idx2) {
        const sensor_pait_t sensor_pair = SensorPairToPairId(sensor_id1, sensor_id2);
        if (SwapSensorPair(sensor_id1, sensor_id2)) {
            pairs_[sensor_pair].insert(index_pair_t(idx2, idx1));
        } else {
            pairs_[sensor_pair].insert(index_pair_t(idx1, idx2));
        }
    }

    void RemoveIndexPair(const sensor_t sensor_id1, const sensor_t sensor_id2, const size_t idx1, const size_t idx2) {
        const sensor_pait_t sensor_pair = SensorPairToPairId(sensor_id1, sensor_id2);
        if (SwapSensorPair(sensor_id1, sensor_id2)) {
            pairs_[sensor_pair].erase(index_pair_t(idx2, idx1));
        } else {
            pairs_[sensor_pair].erase(index_pair_t(idx1, idx2));
        }
    }

    std::vector<index_pair_t> Pairs(const sensor_t sensor_id1, const sensor_t sensor_id2) const {
        const sensor_pait_t sensor_pair = SensorPairToPairId(sensor_id1, sensor_id2);
        const std::unordered_set<index_pair_t>& pairs_set = pairs_.at(sensor_pair);

        std::vector<index_pair_t> pairs(pairs_set.begin(), pairs_set.end());
        if (SwapSensorPair(sensor_id1, sensor_id2)) {
            std::for_each(pairs.begin(), pairs.end(), [](index_pair_t& p){ std::swap(p.idx1, p.idx2); });
        }

        return pairs;
    }

    const std::unordered_map<sensor_pait_t, std::unordered_set<index_pair_t>>& AllPairs() const { return pairs_; }

private:

    std::unordered_map<sensor_pait_t, std::unordered_set<index_pair_t>> pairs_;
};

struct BatchCalibrationOptions {
    // Loss function types: Trivial (non-robust) and Cauchy (robust) loss.
    enum class LossFunctionType { TRIVIAL, CAUCHY };
    LossFunctionType loss_function_type = LossFunctionType::TRIVIAL;

    // Scaling factor determines residual at which robustification takes place.
    double loss_function_scale = 1.0;

    // Wheter to use the additional calibration constraints
    bool use_additional_constraints = true;

    // Whether to print a final summary.
    bool print_summary = true;

    // Ceres-Solver options.
    ceres::Solver::Options solver_options;

    BatchCalibrationOptions() {
      solver_options.function_tolerance = 1e-6;
      solver_options.gradient_tolerance = 1e-10;
      solver_options.parameter_tolerance = 1e-8;
      solver_options.minimizer_progress_to_stdout = false;
      solver_options.max_num_iterations = 100;
      solver_options.max_linear_solver_iterations = 200;
      solver_options.max_num_consecutive_invalid_steps = 5;
      solver_options.max_consecutive_nonmonotonic_steps = 5;
      solver_options.num_threads = -1;
      solver_options.num_linear_solver_threads = -1;
    }

    // Create a new loss function based on the specified options. The caller
    // takes ownership of the loss function.
    ceres::LossFunction* CreateLossFunction() const {
        ceres::LossFunction* loss_function = nullptr;
        switch (loss_function_type) {
        case LossFunctionType::TRIVIAL:
            loss_function = new ceres::TrivialLoss();
            break;
        case LossFunctionType::CAUCHY:
            loss_function = new ceres::CauchyLoss(loss_function_scale);
            break;
        }
        return loss_function;
    }

    inline bool Check() const {
        CHECK_GE(loss_function_scale, 0);
        return true;
    }
};

// Configuration container to setup batch calibration problems
class BatchCalibrationConfig {
public:

    BatchCalibrationConfig()
        : reference_sensor_id_(kInvalidSensorId) { }

    // Access reference sensor id
    inline sensor_t ReferenceSensor() const { return reference_sensor_id_; }
    inline void SetReferenceSensor(const sensor_t sensor_id) { reference_sensor_id_ = sensor_id; }

    // Add / remove sensors from the configuration
    inline void AddSensor(const class Sensor& sensor) {
        const sensor_t sensor_id = sensor.SensorId();
        CHECK(!HasSensor(sensor_id));
        sensors_.emplace(sensor_id, sensor);
    }
    inline bool HasSensor(const sensor_t sensor_id) const { return sensors_.find(sensor_id) != sensors_.end(); }
    inline void RemoveSensor(const sensor_t sensor_id) {
        CHECK(HasSensor(sensor_id));
        sensors_.erase(sensor_id);
    }

    // Access individual sensor information
    inline const class Sensor& Sensor(const sensor_t sensor_id) const { return sensors_.at(sensor_id); }
    inline class Sensor& Sensor(const sensor_t sensor_id) { return sensors_.at(sensor_id); }

    // Access configuration data
    inline const EIGEN_STL_UMAP(sensor_t, class Sensor) & Sensors() const { return sensors_; }
    inline class ObservationDatabase& ObservationDatabase() { return observation_database_; }
    inline std::unordered_set<index_pair_t>& AdditionalConstraints() { return additional_constraints_; }

    inline bool Check() const {
        CHECK(HasSensor(reference_sensor_id_));
        CHECK_NE(reference_sensor_id_, kInvalidSensorId);
        return true;
    }
private:

    sensor_t reference_sensor_id_;

    // sensor_id -> class Sensor
    EIGEN_STL_UMAP(sensor_t, class Sensor) sensors_;

    class ObservationDatabase observation_database_;

    std::unordered_set<index_pair_t> additional_constraints_;

};

inline void PrintSolverSummary(const ceres::Solver::Summary& summary) {
    std::cout << summary.BriefReport() << std::endl;
}

// Nonlinear joint optimization solver.
class BatchCalibration {
public:

    BatchCalibration(const BatchCalibrationOptions& options, const BatchCalibrationConfig& config)
        : options_(options), config_(config) {
        CHECK(options_.Check());
        CHECK(config_.Check());
    }

    inline bool Solve() {
        ceres::LossFunction* loss_function = options_.CreateLossFunction();
        ceres::LocalParameterization* angle_parametrization = AngleLocalParameterization::Create();

        ceres::Problem::Options problem_options;
        ceres::Problem problem(problem_options);

        // Add calibration residuals
        sensor_t reference_id = config_.ReferenceSensor();
        const Sensor& reference_sensor = config_.Sensor(reference_id);
        for (const auto& entry : config_.Sensors()) {
            const sensor_t sensor_id = entry.first;
            Sensor& sensor = config_.Sensor(sensor_id);

            if (sensor_id == reference_id) continue;

            for (const index_pair_t& idx_pair : config_.ObservationDatabase().Pairs(reference_id, sensor_id)) {
                ceres::CostFunction* cost_function = CalibrationResidual::Create(reference_sensor.Obsevation(idx_pair.idx1),
                                                                                 sensor.Obsevation(idx_pair.idx2));
                problem.AddResidualBlock(cost_function, loss_function,
                                         &sensor.Calibration(0), &sensor.Calibration(1), &sensor.Calibration(2), &sensor.Calibration(3));
            }

            // Set angle parameterization
            problem.SetParameterization(&sensor.Calibration(2), angle_parametrization);
        }

        if (options_.use_additional_constraints) {
            // Set additional resuduals
            for (const index_pair_t& sensor_pair : config_.AdditionalConstraints()) {
                const sensor_t sensor1_id = sensor_pair.idx1;
                Sensor& sensor1 = config_.Sensor(sensor1_id);
                const sensor_t sensor2_id = sensor_pair.idx2;
                Sensor& sensor2 = config_.Sensor(sensor2_id);

                for (const index_pair_t& idx_pair : config_.ObservationDatabase().Pairs(sensor1_id, sensor2_id)) {
                    ceres::CostFunction* cost_function = AdditionalCalibrationResidual::Create(sensor1.Obsevation(idx_pair.idx1),
                                                                                               sensor2.Obsevation(idx_pair.idx2));
                    problem.AddResidualBlock(cost_function, loss_function,
                                             &sensor1.Calibration(0), &sensor1.Calibration(1), &sensor1.Calibration(2), &sensor1.Calibration(3),
                                             &sensor2.Calibration(0), &sensor2.Calibration(1), &sensor2.Calibration(2), &sensor2.Calibration(3));
                }
            }
        }

        if (problem.NumResiduals() == 0) {
            return false;
        }

        ceres::Solver::Options solver_options = options_.solver_options;
        solver_options.linear_solver_type = ceres::DENSE_QR;
        solver_options.num_threads = 1;
        solver_options.num_linear_solver_threads = 1;

        std::string solver_error;
        CHECK(solver_options.IsValid(&solver_error)) << solver_error;

        ceres::Solve(solver_options, &problem, &summary_);

        if (solver_options.minimizer_progress_to_stdout) {
            std::cout << std::endl;
        }

        if (options_.print_summary) {
            // TODO: print header
            PrintSolverSummary(summary_);
        }

        return summary_.termination_type == ceres::TerminationType::CONVERGENCE;
    }

    inline const Eigen::Vector4d& Paramerameters(const sensor_t sensor_id) {
        return config_.Sensor(sensor_id).Calibration();
    }

    // Get the Ceres solver summary for the last call to 'Solve'
    inline const ceres::Solver::Summary& Summary() const { return summary_; }

private:

    const BatchCalibrationOptions options_;
    BatchCalibrationConfig config_;
    ceres::Solver::Summary summary_;
};
