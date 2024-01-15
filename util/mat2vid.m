function [Vwrite] = mat2vid(matvid, vidpath, format, framerate, cmap)
%% mat2vid: converts 3D video data in matrix to a video file (.avi, .mp4)
%   Save the video file in the specified directory
%
%   INPUT:
%       matvid   	: matrix contaning video data
%       vidpath     : full file path and name for output video without file extension
%       format      : 'avi' or 'mp4'
%       framerate   : output video frame rate [hz]
%
%   OUTPUT:
%       Vwrite      : video writer object
%

if nargin < 5
    cmap = [];
end

fullpath = [char(vidpath) '.' char(format)];

if strcmp(format, 'avi')
    Vwrite = VideoWriter(fullpath, 'Indexed AVI');
elseif strcmp(format, 'mp4')
    Vwrite = VideoWriter(fullpath, 'MPEG-4');
    Vwrite.Quality = 100;
else
    error('Format must be "avi" or "mp4"')    
end
Vwrite.FrameRate = framerate;

if ~isempty(cmap)
    Vwrite.Colormap = cmap;
end

matvid = squeeze(matvid);
dim = size(matvid); % video size
ch = length(dim); % # of channels in matrix

open(Vwrite)
if ch == 1
    error('Video must be at least 2D')
elseif ch == 2
    warning('Video is 2D: only one frame will be written')
    frame = matvid;
    writeVideo(Vwrite, frame);
elseif ch == 3
    nframe = dim(3);
    for f = 1:nframe
        frame = matvid(:,:,f);
        writeVideo(Vwrite, frame);
    end
elseif ch == 4
    nframe = dim(4);
    for f = 1:nframe
        frame = matvid(:,:,:,f);
        writeVideo(Vwrite, frame);
    end
else
   error('vide must have between 2-4 dimensions') 
end
close(Vwrite)

end