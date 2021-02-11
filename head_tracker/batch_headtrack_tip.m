function [] = batch_headtrack_tip(root, npts, set_mask, loop_mask, playback)
%% batch_headtrack: runs head tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       npts        :   # of points for tip tracker
%       set_mask    :   mask structure or mask mode
%       loop_mask   :   if true, then only set mask for 1st file
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%
%   OUTPUT:
%       -
%

% root = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250\registered';
% npts = 100;
% loop_mask = false;
% set_mask = [1 0];
% playback = 0;

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

headdir = fullfile(PATH,'tracked_head_tip');
mkdir(headdir)
head_mask = set_mask;
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')
    %load(fullfile(PATH,FILES(file)),'regvid','t_v')

    [head_data, head_mask] = track_head(vidData, head_mask, npts, playback);
    
 	save(fullfile(headdir,FILES{file}),'-v7.3','head_data', 'head_mask')
    
    if loop_mask
        % use last mask for next file
    else
        head_mask = set_mask; % make new mask for netx file
    end
                                                                       
end
disp('ALL DONE')
end