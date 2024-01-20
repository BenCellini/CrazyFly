# Head roll tracker guide

The head roll tracker is designed for 2D top or bottom view videos of rigidly tethered flies or videos that have been registered such that the body is fixed in place (see [registration_guide.md](../registration/registration_guide.md)).


## How it works

The [head_roll_tracker.m](head_tracker.m) class will take in a video (must be a MATLAB variable containing a video matrix for now), stabilize the head in the frame, and then estimate the roll angle by comparing the relative widths of the eyes to compute the *roll index*. The roll index can then be calibrated to map to the actual head roll. This technique was originally developed here

    Kim, A. J., Fenk, L. M., Lyu, C. & Maimon, G. Quantitative Predictions Orchestrate Visual Signaling in Drosophila. Cell 168, 280-294.e12 (2017).

where Figure S2 provides more detailed explanation behind the method and calibration.

I implemented this method independently myself (all original code) here

    Cellini, B., Salem, W. & Mongeau, J.-M. Mechanisms of punctuated vision in fly flight. Curr. Biol. 31, 4009-4024.e3 (2021).

Please cite both of these papers if using this specific CrazyFly method.

## Video quality

 Note that for this tracker to work correctly, there must be clear contrast between the eyes and the head:

![head_roll.png](..%2Fimg%2Fhead_roll.png)


## Running
The main class to register videos is [head_roll_tracker.m](head_roll_tracker.m)

    obj = head_roll_tracker(vid, roll_cal)

### Input guide

* `vid` must be a MATLAB variable containing the video as a matrix. This may change in a future release.


* `roll_cal`: roll calibration coefficient. Normally set to 36.33 from Kim et al 2017, but may be best to do your own calibration.

### Output guide

[head_roll_tracker.m](head_roll_tracker.m) will not automatically output a data file like many of the other CrazyFly functions. Instead, the tracking data (and lots of other data) is stored in the `obj` object. Some useful data, like `obj.yaw` & `obj.roll` `obj.roll_idx`, can be saved independent if desired. 

By running

    play_tracking(self, playback, [], savepath, montage_name)

a montage of th roll tracking  can be output.

Example output montage: [example_head_roll_montage.mp4](..%2Fexample_videos%2Ftracked_head%2Fexample_head_roll_montage.mp4)

### Example usage

#### From video file, debugging the initial heading, playing back every frame, & saving the tracking video in the default location

    load('example_videos\example_head_roll.mat', 'vid'); % load .mat video file
    roll_calibration = 36.33; % see Kim et al 2017: "Quantitative Predictions Orchestrate Visual Signaling in Drosophila"
    self = head_roll_tracker(vid, roll_calibration); % run tracker
    
    % Make video montage
    playback = 1;
    savepath = 'example_videos\tracked_head';
    montage_name = 'example_head_roll_montage.mp4';
    play_tracking(self, playback, [], savepath, montage_name)