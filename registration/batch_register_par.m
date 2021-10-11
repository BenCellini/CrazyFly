function [] = batch_register_par(root, vidvar, reg_par, batch_par, crop_xy)
%% batch_register: register selected videos
%
%   INPUT:
%       root        :   root directory
%       vidvar  	:   variable name of video matrix
%       reg_par  	:   use an inner parfor loop for registration if true
%       batch_par  	:   use an outer parfor loop for batch processing if true
%       crop_xy  	:   crop rectangle or boolean to manually crop
%
%   OUTPUT:
%       -
%

if nargin < 5
   crop_xy = []; 
end

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
n_file = length(FILES);

regdir = fullfile(PATH,'registered');
mkdir(regdir)

if ~isempty(crop_xy) && crop_xy
    vid = load(fullfile(PATH,char(FILES(1))), vidvar);
    [~, crop_xy] = imcrop(vid.(vidvar)(:,:,:,1));
    close
end

tic
if batch_par
    parfor (n = 1:n_file, 8)
        register(FILES, PATH, n, regdir, vidvar, reg_par, crop_xy)
    end
else
    for n = 1:n_file
        register(FILES, PATH, n, regdir, vidvar, reg_par, crop_xy)
    end 
end
disp('ALL DONE')
toc
end

%% Single file registration
function [] = register(FILES, PATH, n, regdir, vidvar, reg_par, crop_xy)
    %try
        disp(FILES(n))
        disp('---------------------------------------')
        savepath = fullfile(regdir, FILES{n});
        vidpath = fullfile(PATH, char(FILES(n)));
        vid = load(vidpath, vidvar);
        [regvid, trf] = register_video_par(vid.(vidvar), reg_par, crop_xy);
        save_vid(savepath, regvid, trf, crop_xy)
    %catch
       %warning(['Error registering ' char(FILES(n))])
    %end
end

%% Save for parfor
function [] = save_vid(fpath, regvid, trf, crop_xy)
    save(fpath, '-v7.3', 'regvid', 'trf', 'crop_xy')
end
