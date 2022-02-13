function [] = batch_set_angle(root, vidvar)
%% batch_set_angle: manually set the rotation angle of a video
%
%   INPUT:
%       root     	:   root directory
%       vidvar      :   video variable name
%
%   OUTPUT:
%       -
%

% Select files
[files, path] = uigetfile('*.mat','Select file to switch', root, 'MultiSelect','on');
files = string(files);
n_file = length(files);

% Rotate videos
isWritable = true;
for n = 1:n_file
    disp(n)
    fname = fullfile(path,files(n));
    matObj = matfile(fname, 'Writable', isWritable); % load .mat object
    matObj.rotated = true;
    
    % Manually set rotation angle
    vid = matObj.(vidvar);
    rotangle = set_angle(vid(:,:,1), 'r');
    
    % Rotate video
    [rotvid] = rotate_vid(vid, rotangle);
    matObj.(vidvar) = rotvid;
end
disp('Done')

end