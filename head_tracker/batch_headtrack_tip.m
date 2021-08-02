function [] = batch_headtrack_tip(root, vidvar, npts, mask_mode, loop_mask, playback, export)
%% batch_headtrack: runs head tracker for user selected video files
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
%                       Can also specifiy a mask structure saved from previous head tracker.
%       loop_mask   :   if true, then only set mask for 1st file
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%
%   OUTPUT:
%       -
%

% root = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250\registered';
% vidvar = 'regvid';
% npts = 100;
% mask_mode = [1 0];
% loop_mask = false;
% playback = 0;

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

headdir = fullfile(PATH,'tracked_head_tip');
mkdir(headdir)
head_mask = mask_mode;
for n = 1:nfile
    disp(FILES(n))
    disp('---------------------------------------')
    
    % Load data
    Data = load(fullfile(PATH,FILES(n)));
    t_v = Data.t_v;
    vid = Data.(vidvar);

    if export
    [~,basename] = fileparts(FILES(n));
       vidpath = fullfile(headdir, [char(basename) '.mp4']);
    else
        vidpath = [];
    end
    
    % Run tracker
    [head_data, head_mask] = track_head_vid(vid, head_mask, npts, playback, vidpath);
    
 	save(fullfile(headdir,FILES{n}),'-v7.3','head_data', 'head_mask', 't_v')
    
    if loop_mask
        % use last mask for next file
    else
        head_mask = mask_mode; % make new mask for next file
    end
                                                                       
end
disp('ALL DONE')
end