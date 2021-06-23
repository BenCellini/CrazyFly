function [] = batch_any_tip(root, vidvar, npts, mask_mode, loop_mask, offsetAngle, playback, export)
%% batch_any_tip: runs tip tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       vidvar      :   video variable name
%       npts        :   # of points for tip tracker
%       mask_mode 	:   [x y]. x is a boolean that sets whether initial mask position is 
%                       automatically computed (true) or not (false). y is a boolean that sets 
%                       whether initial mask position is editable by user after being set by 
%                       x (true) or left using the initial position (false). 
%                       Ex: [1 1] (set automatically & and let user edit after)
%                       Can also specifiy a mask structure saved from previous tracker.
%       loop_mask   :   if true, then only set mask for 1st file
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%       export      : (boolean) export video as .mp4
%
%   OUTPUT:
%       -
%

% root = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250_body_fixed';
% vidvar = 'vidData';
% npts = 250;
% mask_mode = [1 1];
% loop_mask = true;
% offsetAngle = [];
% playback = 20;

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
n_file = length(FILES);

bodypart_dir = fullfile(PATH,'tracked_bodypart_tip');
mkdir(bodypart_dir)
bodypart_mask = mask_mode;
for n = 1:n_file
    disp(FILES(n))
    disp('---------------------------------------')
    
    % Load data
    Data = load(fullfile(PATH,FILES(n)));
    %t_v = Data.t_v;
    vid = Data.(vidvar);
    vid = rot90(vid,2);
    [~,basename] = fileparts(FILES(n));
    
    if export
       vidpath = fullfile(bodypart_dir, [char(basename) '.mp4']);
    else
        vidpath = [];
    end
    
    % Run tracker
    [bodypart, bodypart_mask, offsetAngle] = track_any(vid, bodypart_mask, ...
        npts, offsetAngle, playback, vidpath);
    
 	save(fullfile(bodypart_dir,FILES{n}),'-v7.3','bodypart', 'bodypart_mask', 'offsetAngle')
    
    if loop_mask
        % use last mask for next file
    else
        bodypart_mask = mask_mode; % make new mask for netx file
    end
                                                                       
end
disp('ALL DONE')
end