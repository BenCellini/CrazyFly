# CrazyFly

## Overview

CrazyFly is a set of tools written in MATLAB for pulling out kinematics from videos of rigidly and magnetically tethered  (fly body can rotate about yaw axis) fruit flies. While not tested thoroughly on videos of other organisms, CrazyFly should generalize fairly well.

Please cite:

    Cellini, B., Salem, W. & Mongeau, J.M. Complementary feedback control enables effective gaze stabilization in animals. Proc. Natl. Acad. Sci. 119, e2121660119 (2022).

and link this repository if using CrazyFly software in published works.

The latest version of CrazyFly has been tested on MATLAB 2022a, but should work for most newer releases. Some older releases have slightly different syntax for handling movable points on figures that might throw some errors.

CrazyFly should work really well for videos of similar quality to the example videos provided in [example_videos](example_videos). I have analyzed thousands of videos with this software suite and my best estimate of the success rate (correctly tracked kinematics) would be 99%. There are of course possible bugs and other things I did not account for when analyzing videos that deviate from the quality of my videos, so don't expect perfect results if your videos are not a close match to mine.

Any questions about CrazyFly or suggestions for new features/improvements can be directed to me (Benjamin Cellini): `bcellini00@gmail.com`. Please include "CrazyFly" in the subject line. I am happy to help troubleshoot when I am not too busy.

## Tips on getting good videos
All the CrazyFly functions work infinitely better with high quality videos, and it is worth the extra time spent optimizing video quality before attempting any kinematic tracking. Here are some tips for taking videos of tethered fruit flies:

### Frame rate

A frame rate around 100 fps is ideal if you plan on measuring wing stroke angles, as each image will contain approximately 2 wing strokes, resulting in a nice blurred wing with a defined edge. This makes extracting the wing envelope edge a much easier task than with higher frame rates. Higher frame rates may be beneficial for tracking the head yaw angle, but 100 fps is still generally enough.

### Camera angle

Choose an angle that is aligned with the wing stroke plane and head yaw plane as possible. The wing stroke edge should appear high and close to the fly head, and you should be able to see a clear gap between the fly head and neck. The base of the antennae on the head should also be visible. This is critical, as CrazyFly relies on the antennae features to measure the head yaw angle.

### Lighting
CrazyFly is designed to work wih front lighting, meaning that the fly is the only illuminated object (appears white) in the frame and everything else is (ideally) black. Backlit videos should also work, but have to be inverted. CrazyFly does this automatically, but this has not been thoroughly tested. In either case, having enough lighting + contrast to clearly resolve edges is very helpful.

## Example videos

Examples of good quality videos are provided files in the [example_videos](example_videos) directory as both `.mp4` and `.mat` files. Most CrazyFly functions can take in either path's to video files (`.mp4`, `.avi`, etc.) or videos stored in a matrix data (from `.mat` files). If your videos are stored as `.mat` files and are too large to be loaded into MATLAB, I recommend converting your files to `.mp4` before using any CrazyFly functions.
The [mat2vid.m](util/mat2vid.m) &  [vid2mat.m](util/vid2mat.m) functions I provide may be useful for converting your files between formats.

## Main CrazyFly functions
Each function comes with its own guide file. Make sure to have the path to the main tracker function (e.g., [body_tracker.m](body_tracker/body_tracker.m)) and the [util](util) directory  on the MATLAB path before running any functions. The file [run_trackers.m](run_trackers.m) shows example usage for all main CrazyFly functions, but it is recommended to read the guides for all functions (linked below) first.

### Registration

 [registration_guide.md](registration/registration_guide.md)

 Register videos of magnetically tethered flies such that every frame is shifted to make it appear that the body does not move. This is a useful preprocessing step that makes further analysis (for example, head or abdomen tracking) much simpler. Skip this if you are working with videos of rigidly tethered (body-fixed) flies.

### Body yaw tracking

[body_tracker_guide.md](body_tracker/body_tracker_guide.md)

Tracks the body yaw angle of magnetically tethered flies. Note that registering your videos also gives you the body angle, but this method is much faster. Skip this if you are working with videos of rigidly tethered flies.

### Head yaw tracking

[head_tracker_guide.md](head_tracker/head_tracker_guide.md)

Tracks the head yaw angle with respect to the neck joint of rigidly tethered flies (or registered videos of magnetically tethered flies).

### Head roll tracking

[head_roll_tracker_guide.md](head_tracker/head_roll_tracker_guide.md)

Attempts to estimate head roll by from a single top or bottom 2D camera view by computing the ratio of eye widths. Requires higher contrast videos like [example_head_roll.mp4](example_videos/example_head_roll.mp4).


### Abdomen yaw tracking

[abdomen_tracker_guide.md](abdomen_tracker/abdomen_tracker_guide.md)

Tracks the abdomen yaw angle with respect to the abdomen joint of rigidly tethered flies (or registered videos of magnetically tethered flies).

### Wing beat amplitude tracking

Not yet supported. I recommend https://github.com/BenCellini/Benifly for wing tracking.

## Notes

* If your videos are stored in `.mat` files and are too big to load into ram, convert them to `.mp4` files first. [mat2vid.m](util/mat2vid.m) is a useful function for this purpose.

## Disclaimer
*Note that CrazyFly has nothing to do with the Drone platform Crazyflie (https://www.bitcraze.io/products/crazyflie-2-1/), I just didn't realize there was anything else called Crazyfly when I made this repository...but it's too late to change the name.*