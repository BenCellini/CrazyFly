function [vid] = vid2mat(vidpath, savedir)
%% vid2mat: converts video file (*.mp4, *.avi, etc.) to matlab matrix
%
%   INPUT:
%       vidpath         : full path to video file
%       savedir         : directory to save video to, don't save if empty
%
%   OUTPUT:
%       vid             : video matrix
%

if nargin < 2
    savedir = [];
end

% Create video-reader object
vidpath = char(vidpath);
Vread = VideoReader(vidpath);

% Preallocate video matrix
init_frame = read(Vread, 1);
dim = size(init_frame);
vid = uint8(zeros([dim(1:2) Vread.NumFrames]));

% Read in and store all frames
for n = 1:Vread.NumFrames
    frame = read(Vread, n);
    frame_gray = uint8(im2gray(frame));
    vid(:, :, n) = frame_gray;
end

% Save .mat file containing video if speciifed
if ~isempty(savedir)
    [~,basename,~] = fileparts(vidpath);
    matpath = fullfile(savedir, [basename '.mat']);
    save(matpath, 'vid', '-v7.3')
    disp(['Saved: ' matpath])
end

end