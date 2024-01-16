function [viddata] = get_vid_data(vid, save_path, mat_var_name, file_add)
%% get_vid_data: parse the 'vid' variable
%
%   INPUT:
%       vid             : raw grayscale video data in matrix form
%                         or a path to a video file (.mp4 or .avi)
%       save_path       : path to save new data
%                         default is 'data' folder in 'vid' directory
%       mat_var_name    : name of the matlab variable containing the video data
%                         only required if inputing vid_path to .mat file
%       file_add        : string to add at end of output file name
%                         also name of new folder if using defualt save location
%
%   OUTPUT
%       viddata.
%           matflag         : (boolean) if 'vid' is a matrix or a path to a video file
%           first_frame     : first frame from video
%           n_frame         :
%           Vread           : video-reader object
%                             if 'vid' is a matrix or a path to a .mat file, then returns empty
%           vid             : video matrix
%                             if 'vid' is a path to a video file, then returns empty
%

% Parse inputs
if nargin < 4
    file_add = '';
    if nargin < 3
        mat_var_name = 'vid';
        if nargin < 2
            save_path = [];
        end
    end
end

% Check how video was input
if ischar(vid) || isstring(vid)  % video path was input
    [vid_dir, vid_name, vid_ext] = fileparts(vid); % parse video path

    % Check if .mat file was given, otherwise assume video file format (.mp4, .avi, etc)
    if strcmp(vid_ext, '.mat')
        matflag = true; % video stored in MATLAB variable
        vid_data = load(vid, mat_var_name); % load video variable from .mat file
        vid = vid_data.(mat_var_name); % store video data
        n_frame = size(vid, 3); % # of frames
        first_frame = vid(:, :, 1); % get 1st frame
    else
        matflag = false; % video NOT stored in MATLAB variable, need to read from file instead
        Vread = VideoReader(vid); % video-reader object
        n_frame = Vread.NumFrames; % # of frames
        first_frame = im2gray(read(Vread, 1)); % get 1st frame
    end

else % video matrix was input
    matflag = true;
    vid = squeeze(vid); % in case it's 4D
    n_frame = size(vid, 3); % # of frames
    first_frame = vid(:, :, 1); % get 1st frame

    if isempty(save_path)
        error('Must specify path to save registered video when inputing video as matrix.')
    end
end

% Set path to save registered video, always save as .mp4
if isempty(save_path)
    vid_dir = fullfile(vid_dir, file_add);
    save_path = fullfile(vid_dir, [char(vid_name) '_' file_add '.mp4']);
    save_data_path = fullfile(vid_dir, [char(vid_name) '_' file_add '_data.mat']);
else
    [vid_dir, vid_name, ~] = fileparts(save_path); % parse path
    save_data_path = fullfile(vid_dir, [char(vid_name) '_' file_add '_data.mat']);
end

% Make output folder if it does not already exist
[vid_dir, ~, ~] = fileparts(save_path);
if ~exist(vid_dir, 'dir') && ~isempty(vid_dir)
   mkdir(vid_dir)
end

% Check lighting
pixel_mean = mean(first_frame, 'all');
if pixel_mean > 100 % back lighting detected
    first_frame = imcomplement(first_frame);
    invertflag = true;
else
    invertflag = false;
end

% Store data
viddata.matflag = matflag;
viddata.first_frame = first_frame;
viddata.invertflag = invertflag;
viddata.n_frame = n_frame;

if matflag
    viddata.Vread = [];
    viddata.vid = vid;
else
    viddata.Vread = Vread;
    viddata.vid = [];
end

viddata.save_vid_path = save_path;
viddata.save_data_path = save_data_path;

end
