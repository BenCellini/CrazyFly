function [] = batch_bodytrack(root, playback, heading_debug, par, export)
%% batch_bodytrack: runs body tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       playback    :   playback rate (show a frame in increments of "playback").
%                       If false, then onlt show the 1st frame. (default = 1)
%       head_debug  :   always bring up the heading angle check window if true
%       par         :   use parallel processing
%
%   OUTPUT:
%       -
%

% playback = 10;
% root = 'H:\EXPERIMENTS\MAGNO\Experiment_SOS';

if nargin < 5
    export = false;
    if nargin < 4
        par = false; % default
        if nargin < 3
            heading_debug = false; % default
            if nargin < 2
                playback = 1; % default
                if ~nargin
                    root = ''; % root is current folder
                end
            end
        end
    end
end

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

bodydir = fullfile(PATH,'tracked_body');
mkdir(bodydir)
for n = 1:nfile
    disp(FILES(n))
    disp('---------------------------------------')
    load(fullfile(PATH,char(FILES(n))),'vidData','t_v')
    [~,basename] = fileparts(FILES(n));
   	close all
    
    if export
       vidpath = fullfile(bodydir, [char(basename) '.mp4']);
    else
        vidpath = [];
    end
    
    [bAngles,imgstats,initframe] = bodytracker_vid(vidData, playback, heading_debug, par, vidpath);

    save(fullfile(bodydir,FILES{n}),'-v7.3','bAngles','imgstats','initframe','t_v')
end
disp('ALL DONE')
end