function [] = batch_headtrack_roll(root_yaw)
%% batch_headtrack_roll: runs head roll tracker for user selected video files
%
%   INPUT:
%       root_yaw   	:   root directory with yaw angles
%       root_vid   	:   root directory with video angles

%       npts        :   # of points for tip tracker
%       set_mask    :   mask structure or mask mode
%       loop_mask   :   if true, then only set mask for 1st file
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%
%   OUTPUT:
%       -
%

[FILES, path_yaw] = uigetfile({'*.mat', 'MAT-files'},'Select yaw angles', root_yaw, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

sub_fold = regexp(root_yaw, '\', 'split');
path_vid = fullfile(sub_fold{1:end-1});
%path_vid = root_vid;b

eye_ratio = 0.3;
roll_cal = 36.33;

headdir = fullfile(path_yaw,'tracked_head_roll');
mkdir(headdir)
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(path_yaw,FILES(file)),'head_data','head_mask')
    load(fullfile(path_vid,FILES(file)),'regvid','t_v')

    [roll, roll_idx] = track_head_roll(regvid, head_data.angle, head_mask.move_points.rot, ...
        eye_ratio, roll_cal, true);
    
 	save(fullfile(headdir,FILES{file}),'-v7.3','roll', 'roll_idx', 't_v')
                                                                       
end
disp('ALL DONE')
end