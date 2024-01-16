# CrazyFly

## Overview

CrazyFly is a set of tools written in MATLAB for pulling out kinematics from videos of rigidly and magnetically tethered  (fly body can rotate about yaw axis) fruit flies. While not tested thoroughly, CrazyFly can also work well with other insects (such as tethered hawk moths).

Please cite:

    Cellini, B., Salem, W. & Mongeau, J.M. Complementary feedback control enables effective gaze stabilization in animals. Proc. Natl. Acad. Sci. 119, e2121660119 (2022).

and link this repository if using CrazyFly software in published works.

Any questions about CrazyFly can be directed to me (Benjamin Cellini): `bcellini00@gmail.com`.

## Tips on getting good videos
All the CrazyFly functions work infinity better with high quality videos and it is worth the extra time spent optimizing video quality before attempting any kinematic tracking. Here are some tips for taking videos of tethered fruit flies:

### Frame rate

A frame rate around 100 fps is ideal if you plan on measuring wing stroke angles, as each image will contain approximately 2 wing strokes, resulting in a nice blurred wing with a defined edge. This makes extracting the wing envelope edge a much easier task than with higher frame rates. Higher frame rates may be beneficial for tracking the head yaw angle, but 100 fps is still enough to capture even faster head movements.

### Camera angle

Choose an angle that is as in line with the wing stroke plane and head yaw plane as possible. The wing stroke edge should appear high and close to the fly head and you should be able to see a clear gap between the fly head and neck. The base of the antennae on the head should also be visible. This is critical, as CrazyFly relies on the antennae features to measure the hed yaw angle.

### Lighting
CrazyFly is designed to work wih front lighting, meaning that the fly is the only illuminated object (appears white) in the frame and everything else is (ideally) black. Backlit videos should also work, but have to be inverted. CrazyFly does this automatically, but this has not been thoroughly tested. In either case, having enough lighting + contrast to clearly resolve edges is very helpful.


## Example videos

Examples of good quality videos are provided files in the [example_vidoes](example_vidoes) directory as both `.mp4` and `.mat` files. Most CrazyFly functions can take in either path's to video files (`.mp4`, `.avi`, etc.) or videos stored in a matrix data (from `.mat` files). If your videos are stored as `.mat` files and are too large to be loaded into MATLAB, I recommend converting your files to `.mp4` before using any CrazyFly functions.
The [mat2vid.m](util%2Fmat2vid.m) &  [vid2mat.m](util%2Fvid2mat.m) functions I provide may be useful for converting your files between formats.

## Main CrazyFly functions
Each function comes with its own guide file. Make sure to have the path to the main tracker function (e.g., [body_tracker.m](body_tracker%2Fbody_tracker.m)) and the [tools](tools) directory  on the MATLAB path before running any functions. The file [run_trackers.m](run_trackers.m) shows example usage for all main CrazyFly functions, but it is recommended to read the guides for all functions (linked below) first.

### Registration

 [registration_guide.md](registration%2Fregistration_guide.md)

 Register videos of magnetically tethered flies such that every frame is shifted to make it appear that the body does not move. This is a useful preprocessing step that makes further analysis (for example, head or abdomen tracking) much simpler. Skip this if you are working with videos of rigidly tethered (body-fixed) flies.

### Body yaw tracking

[body_tracker_guide.md](body_tracker%2Fbody_tracker_guide.md)

Tracks the body yaw angle of magnetically tethered flies. Note that registering your videos also gives you the body angle, but this method is much faster. Skip this if you are working with videos of rigidly tethered flies.

### Head yaw tracking

[head_tracker_guide.md](head_tracker%2Fhead_tracker_guide.md)

Tracks the head yaw angle with respect to the neck joint of rigidly tethered flies (or registered videos of magnetically tethered flies).

### Abdomen yaw tracking

[abdomen_tracker_guide.md](abdomen_tracker%2Fabdomen_tracker_guide.md)

Tracks the abdomen yaw angle with respect to the abdomen joint of rigidly tethered flies (or registered videos of magnetically tethered flies).

### Head roll tracking


## Notes
* Many of the functions 

## Disclaimer
*Note that CrazyFly has nothing to do with the Drone platform Crazyflie (https://www.bitcraze.io/products/crazyflie-2-1/), I just didn't realize there was anything else called Crazyfly when I made this repository...but it's too late to change the name.*