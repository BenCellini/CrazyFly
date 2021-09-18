function [] = batch_register_par(root, crop_xy)
%% batch_register: register selected videos
%
%   INPUT:
%       root        :   root directory
%       crop_xy  	:   crop rectangle or boolean to manually crop
%
%   OUTPUT:
%       -
%

if nargin < 2
   crop_xy = []; 
end

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
n_file = length(FILES);

regdir = fullfile(PATH,'registered');
mkdir(regdir)

if length(crop_xy) == 1
    load(fullfile(PATH,char(FILES(1))), 'video')
    [~, crop_xy] = imcrop(video(:,:,:,1));
    close
end

parfor n = 1:n_file
    try
        disp(FILES(n))
        disp('---------------------------------------')
        data = load(fullfile(PATH,char(FILES(n))), 'video');

        [regvid,trf] = register_video_par(data.video, crop_xy);
        fpath = fullfile(regdir,FILES{n});
        save_vid(fpath, regvid, trf, crop_xy)
    catch
       warning(['Error registering ' char(FILES(n))]) 
    end
end
disp('ALL DONE')
end

%% save fucntion for parfor
function [] = save_vid(fpath, regvid, trf, crop_xy)
    save(fpath,'-v7.3','regvid','trf','crop_xy')
end






