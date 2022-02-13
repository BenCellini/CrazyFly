function [rotvid] = rotate_vid(vid, angles, method, bbox)
%% rotate_vid: rotate every frame in video
%
%   INPUT:
%       vid    	: .mat video matrix
%       angles 	: rotation angle(s)
%
%   OUTPUT:
%       rotvid 	: rotated video
%

if nargin < 4
    bbox = 'crop';
    if nargin < 3
        method = 'nearest';
    end
end

vid = squeeze(vid);
dim = size(vid);

if length(angles) == dim(3) % set angle for each frame
    % do nothing
elseif length(angles) == 1 % same angle for every frame
    angles = repmat(angles, [dim(3) 1]);
else
    error('''angles'' must be a scalar or a vector of length size(vid,3)')
end

% Rotate each frame
for f = 1:dim(3) % each frame
    frame = vid(:,:,f);
    rotframe = imrotate(frame, angles(f), method, bbox); % rotate
    rotvid(:,:,f) = rotframe;
end

end