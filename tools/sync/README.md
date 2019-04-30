# Synchronization tool
Computes synchronized incremental motions from trajectories estimated for multiple sensors.

## 1. Dependencies

See main project [README.md](../../README.md)

## 2. Data Format

The input are estimated 2D trajectories from asynchronous sensors. The estimated trajectory are represented in a TUM-like format file. Each trajecotry is described in a text file, where each line is in the format:
```
timestamp x y yaw
```
where the `timestamp` is in seconds, and the `yaw` angle in radians.

## 3. Usage

The synchronizer can be invoked with the following syntax:
```
sync [options] <trajectory1> <trajectory2> ...
```

Notes:

* The incremental motions are synchronized with respect to a sensor as time reference, defined to be the first input trajectory for convenience.

* The output are synchronized incremental motions, in independent files for each sensor, with the output format required by the [2D calibration implementation](https://github.com/dzunigan/calibration2d). The incremental motions are stored in `output_dir/1.txt`, `output_dir/2.txt`, etc. The default `output_dir` option value will store the synchronized motions in the currect directory.

* Use the `-s n` option to control the output rate. This option sets the *skip* value for the reference trajectory (default to 0). This means that if the reference frame rate is 30 Hz, with `-s 3` the synchronization rate will be about 30 / (3+1) = 7.5 Hz.

## 4. Example

See project [README.md](../../README.md)
