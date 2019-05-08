#! /bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: simulation.sh <noise_level>"
fi

noise_level=$1

lin_std=$(bc <<<"0.001*$noise_level")
ang_std=$(bc <<<"0.03*$noise_level")
point_std=$(bc <<<"0.01*$noise_level")

./build/simulate3d --lin_std=$lin_std --ang_std=$ang_std odo_$noise_level.txt cam_$noise_level.txt
if [ $? -ne 0 ]; then exit $?; fi

./build/simulate_ground --point_std=$point_std points_$noise_level.txt
if [ $? -ne 0 ]; then exit $?; fi

ground_calib=$(./build/ground_calib --csv points_$noise_level.txt)
if [ $? -ne 0 ]; then exit $?; fi

ground_calib=($ground_calib)

./build/tools/planar/planar odo_$noise_level.txt odo_${noise_level}_planar.txt
if [ $? -ne 0 ]; then exit $?; fi

./build/tools/planar/planar --pitch=${ground_calib[1]} --roll=${ground_calib[2]} cam_$noise_level.txt cam_${noise_level}_planar.txt
if [ $? -ne 0 ]; then exit $?; fi

./build/tools/sync/sync odo_${noise_level}_planar.txt cam_${noise_level}_planar.txt
if [ $? -ne 0 ]; then exit $?; fi

motion_calib=$(./build/calibrate 1.txt 2.txt)
if [ $? -ne 0 ]; then exit $?; fi

motion_calib=($motion_calib)

echo "${motion_calib[0]} ${motion_calib[1]} ${ground_calib[0]} ${motion_calib[2]} ${ground_calib[1]} ${ground_calib[2]} ${motion_calib[3]}"

