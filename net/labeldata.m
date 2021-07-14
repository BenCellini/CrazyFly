function [rc] = labeldata(root, target, fpv, vpf)
%% labeldata: get videos and pull out frames labeled as head-up or head-down, store in new directory
%
%   INPUT:
%       root    	:   root directory to load raw images
%       target      :   target directory to save processed images
%       fpv         :   how many frames per video to save
%       vps         :   
%
%   OUTPUT:
%       rc          :   
%

% root = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_amp_3.75';
% root = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v2';
% root = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250';
% root = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_amp_v3';
% target = 'Q:\Research\fly_image';
% fpv = 10;
% vpf = 22;

rng(1) % for reproducability

[FILES,PATH] = uigetfile({'*.mat'},'Select head-free data', root, 'MultiSelect','on');
FILES = cellstr(FILES);
n_file = length(FILES);

% Make directories to store images and flipped images
updir = fullfile(target, 'Up');
downdir = fullfile(target, 'Down');
mkdir(updir)
mkdir(downdir)

% Pick videos to use from all slected files
if ~isempty(vpf)
    rand_vids = randperm(n_file);
    rand_vids = rand_vids(1:vpf);
else
    rand_vids = 1:n_file;
end
n_vid = length(rand_vids);

[~,expname,~] = fileparts([PATH(1:end-1) '.t']);
rc = cell(n_vid,1);
for n = 1:n_vid
    disp(n)
    load(fullfile(PATH, FILES{n}), 'vidData');
    [~,basename,~] = fileparts(FILES{n});
    
    vid = squeeze(vidData);
    dim = size(vid);
    
    % Pick frames to use from video
    if ~isempty(fpv)
        rand_frames = randperm(dim(3));
        rand_frames = rand_frames(1:fpv);
    else
        rand_frames = 1:dim(3);
    end
    n_frame = length(rand_frames);
    
    % Get the fly ROI in each image and try to find the right heading
    imgname = [expname '_' basename];
    rc{n} = zeros(n_frame,2);
    for f = 1:n_frame
        raw_frame = vid(:,:,rand_frames(f));
        [~,~,fly_frame,~] = getflyroi(raw_frame, [221 132], 0.25, 2, 0.12);
        rc{n}(f,:) = size(fly_frame); % row & column pixel size of fly frame
        framename = [imgname '_frame_' num2str(rand_frames(f)) '.jpg'];
        
        flip_frame = rot90(fly_frame,2); % flip frame upside down
        imwrite(fly_frame, fullfile(updir, framename))
        imwrite(flip_frame, fullfile(downdir, framename))
    end
    
end
rc_all = cat(1,rc{:});

end