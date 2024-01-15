function [trf, body_angles_from_trf] = register_video(vid, save_vid_path, mat_var_name, debug_heading, flip_vid)
%% register_video: register each frame of a video to the same reference
%
%   INPUT:
%       vid             : raw grayscale video data in matrix form
%                         or a path to a video file (.mp4 or .avi)
%       reg_vid_path    : path to save registered video
%                         default is 'registered' folder in 'vid' directory
%       mat_var_name    : name of the matlab variable containing the video data
%                         only required if inputing vid_path to .mat file
%       flip_vid        : if true, then flip the registered video from left to right 
%                         (usually only set to true if a bottom view video)
%       debug           : show debug figure for setting initial heading
%                         0=never, 1=always, 2=if close call (default = 0)
%
%   OUTPUT:
%       trf   	              : 2D affine transformation for each frame
%       body_angles_from_trf  : body angles calcualted from registration transforma                         tions in 'trf'
%

% Parse inputs
if nargin < 5
    flip_vid = false;
    if nargin < 4
        debug_heading = 0;
        if nargin < 3
            mat_var_name = 'vid';
            if nargin < 2
                save_vid_path = [];
            end
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
    
    if isempty(save_vid_path)
        error('Must specify path to save registered video when inputing matrix as video.')
    end
end

% Set path to save registered video, always save as .mp4
if isempty(save_vid_path)
    vid_dir = fullfile(vid_dir, 'registered');
    save_vid_path = fullfile(vid_dir, [char(vid_name) '_registered.mp4']);
    save_data_path = fullfile(vid_dir, [char(vid_name) '_registered_data.mat']);
else
    [vid_dir, vid_name, ~] = fileparts(save_vid_path); % parse path
    save_data_path = fullfile(vid_dir, [char(vid_name) '_registered_data.mat']);
end

% Make output folder if it does not already exist
[vid_dir, ~, ~] = fileparts(save_vid_path);
if ~exist(vid_dir, 'dir') && ~isempty(vid_dir)
   mkdir(vid_dir)
end

% Get intial heading orientation
[refangle, flip] = find_heading(first_frame, debug_heading);
if ~flip
    refangle = refangle - 180;
end

% Set up optimizer for registration
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 150;
optimizer.MaximumStepLength = 0.03;
optimizer.MinimumStepLength = 0.0002;

% Registration setup
trf = cell(n_frame,1); % store 2D affine transformations here
fixed = double(first_frame); % initial frame to register to
sz = imref2d(size(fixed));

% Get the initial transform for the 1st frame
trf{1} = imregtform(fixed, fixed, 'rigid', optimizer, metric);

% Create video writer object
Vwrite = VideoWriter(save_vid_path, 'MPEG-4');
Vwrite.Quality = 100;
Vwrite.FrameRate = 60;
open(Vwrite)

% Register each frame with respect to the 1st frame
tic
for n = 1:n_frame
    if mod(n, 10) || (n == 1) || (n== n_frame)
        fprintf([int2str(n) '\n'])
    end

    % Get frame
    if matflag % from matrix
        frame = vid(:, :, n);
    else % from video reader 
        frame = im2gray(read(Vread, n));
    end

    % Convert frame to double
    frame = double(frame);

    % Find the last rigid transform
    for jj = n-1:-1:1
        if(isRigid(trf{jj}))
            break
        end
    end

    % Find new transform
    if n > 1
        trf{n} = imregtform(frame, fixed, 'rigid', optimizer, metric,...
                            'InitialTransformation', trf{jj});
    end

    % Warp current frame by new transform
    reg = imwarp(frame, trf{n}, 'OutputView', sz);

    % Rotate so head is facing up
    reg_rot = uint8(imrotate(reg, 90-refangle, 'crop'));

    % Write regsitered & rotated frame to output video
    if flip_vid
        reg_rot = fliplr(reg_rot);
    end
    writeVideo(Vwrite, reg_rot);

    % Update the fixed frame
    fixed = (fixed*n + reg) / (n + 1);
end

% Calculate body angles from transformation data
T11 = cellfun(@(x) x.T(1,1), trf);
T12 = cellfun(@(x) x.T(1,2), trf);
body_angles_from_trf = rad2deg(unwrap(atan2(T12, T11))) + (90 - refangle);

% Save transformation data & body angles
save(save_data_path, 'trf', 'body_angles_from_trf', '-v7.3')

disp('DONE')
toc
end