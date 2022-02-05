# Automatic Multi-Sensor Extrinsic Calibration for Mobile Robots
Implementation of an automatic method to estimate the extrinsic calibration parameters of multiple sensors, explicitly handling the scale ambiguity arising from monocular cameras.

**Authors:** [David Zuñiga-Noël](http://mapir.uma.es/dzuniga), [Jose-Raul Ruiz-Sarmiento](https://mapir.isa.uma.es/mapirwebsite/?p=1366), [Ruben Gomez-Ojeda](https://mapir.isa.uma.es/mapirwebsite/?p=2155), [Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536)

**Related Publications:** [Automatic Multi-Sensor Extrinsic Calibration for Mobile Robots](https://arxiv.org/abs/1906.04670)

**License:**  [GPLv3](https://raw.githubusercontent.com/dzunigan/robot_autocalibration/master/LICENSE.txt)

## 1. Dependencies

* Boost (1.58.0.1ubuntu1)
   ```
   sudo apt install libboost-all-dev
   ```
* Ceres (1.14.0-facb199)

   See ceres [documentation](http://ceres-solver.org/installation.html#linux).
   
* CMake (3.5.1-1ubuntu1)
   ```
   sudo apt install cmake
   ```
* Eigen (3.3~beta1-2)
   ```
   sudo apt install libeigen3-dev
   ```
* Gflags (2.1.2-3)
   ```
   sudo apt install libgflags-dev
   ```
* Glog (0.3.4-0.1)
   ```
   sudo apt install libgoogle-glog-dev
   ```

## 2. Data Format

The input for the ground plane calibration are 3D point observations from the ground, in local sensor coordinates. The input is a csv text file as:
```
x,y,z
```

The input for the  are synchronized per-sensor incremental 2D (coplanar) motions. For each sensor, incremental motions are grouped into a single text file, where each line contains a synchronized 2D incremental motion (angles in radians) in the following format:
```
x y yaw
```

Note: all sensors must have the same number of observations.

## 3. Usage

The multisensor calibration method can be invoked as:
```
calibrate [options] <motions1> <motions2> ...
```

Notes:

* The method calibrates all sensors with respect to a reference one. The reference sensor is assumed to be the one observing the first incremetal motions.

* Use the `--scale_ambiguous=n,m,...` option to indicate which motions have scale ambiguity.

## 4. Example


