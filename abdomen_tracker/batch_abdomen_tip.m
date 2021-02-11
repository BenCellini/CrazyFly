function [] = batch_abdomen_tip(root, npts, set_mask, loop_mask, playback)
%% batch_abdomen_tip: runs abdomen tracker for user selected video files
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

% root = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250_body_fixed';
% npts = 150;
% loop_mask = false;
% set_mask = [1 1];
% playback = 20;

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

abdomen_dir = fullfile(PATH,'tracked_abdomen_tip');
mkdir(abdomen_dir)
abdomen_mask = set_mask;
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    %load(fullfile(PATH,FILES(file)),'vidData','t_v')
    load(fullfile(PATH,FILES(file)),'regvid','t_v')

    [abdomen_data, abdomen_mask] = track_abdomen(regvid, abdomen_mask, npts, playback);
    
 	save(fullfile(abdomen_dir,FILES{file}),'-v7.3','abdomen_data', 'abdomen_mask')
    
    if loop_mask
        % use last mask for next file
    else
        abdomen_mask = set_mask; % make new mask for netx file
    end
                                                                       
end
disp('ALL DONE')
end