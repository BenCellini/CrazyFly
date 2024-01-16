%% Examples of main CrazyFly functions
%
% Make sure to start in the root of this repository & that
% all subfolders in the repository folder are on the MATLAB path.
%
% WARNING: running these functions may overwrite some of the vidoes & 
% data files already in the 'example_vidoes' diectory.
%

addpath(genpath(pwd)) % add the repo to the path

clear ; close all ; clc

%% Registration
% Register the example video of a magnetically tethered (body-free) fly 
% in the 'example_vidoes' directory. Save to the default path: '..\registered'.
clear
vid = 'example_vidoes\example_body_free.mp4'; % path to video file
save_vid_path = []; % path to new file (e.g., 'registered_video.mp4'). Leaving empty uses defualt path
mat_var_name = []; % only required if loading vidoe from a .mat file
debug_heading = 1; % if 1 then bring up the debug heading window automatically
flip_vid = false; % do not flip the video if false

[trf, body_angles_from_trf] = register_video(vid, save_vid_path, mat_var_name, debug_heading, flip_vid);

plot(body_angles_from_trf)

%% Body tracking
% Track the body in the example video of a magnetically tethered (body-free) fly 
% in the 'example_vidoes' directory. Save to the default path: '..\tracked_body'.
clear
vid = 'example_vidoes\example_body_free.mat'; % path to video file
save_path = []; % path to new file (e.g., 'registered_video.mp4'). Leaving empty uses defualt path.
mat_var_name = []; % only required if loading vidoe from a .mat file
debug_heading = 1; % if 1 then bring up the debug heading window automatically
flip_vid = false; % do not flip the video if false
playback = 1; % playback rate (1 = every frame)
vidpath = []; % path to save trackign video. Leaving empty uses defualt path

[angle, imgstats, initframe] = body_tracker(vid, save_path, mat_var_name, debug_heading, flip_vid, playback, vidpath);

plot(angle)

%% Body tracking
% Track the body in the example video of a magnetically tethered (body-free) fly 
% in the 'example_vidoes' directory. Save to the default path: '..\tracked_body'.
save_path = 'test.mat'; % path to new file (e.g., 'registered_video.mp4'). Leaving empty uses defualt path.
mat_var_name = 'vid'; % only required if loading vidoe from a .mat file
debug_heading = 1; % if 1 then bring up the debug heading window automatically
flip_vid = false; % do not flip the video if false
playback = 1; % playback rate (1 = every frame)
vidpath = false; % path to save trackign video. Leaving empty uses defualt path

[angle, imgstats, initframe] = body_tracker(vid_back, save_path, mat_var_name, debug_heading, flip_vid, playback, vidpath);

plot(angle)

%% Comparison of body angles from registration and body-tracker
clear ; close all

load('example_vidoes\registered\example_body_free_registered_data.mat')
load('example_vidoes\tracked_body\example_body_free_tracked_body_data.mat')

hold on;
plot(body_angles_from_trf, 'k', 'LineWidth', 1)
plot(angle, 'r', 'LineWidth', 1)
legend('from registration', 'from body-tracker')

%% Head tracking
% Track the body in the example registered video of a magnetically tethered (body-free) fly 
% in the 'example_vidoes\registered' directory. Save to the default path: '..\tracked_head'.





