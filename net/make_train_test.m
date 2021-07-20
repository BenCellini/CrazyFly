function [] = make_train_test(root, target, splitp)
%% make_train_test: splits folder of labeled images into train and test data
%
%   INPUT:
%       root    	:   root directory to load raw images
%       target      :   target directory to save processed images
%       splitp   	:   train/test split percent
%
%   OUTPUT:
%       -
%

% root = 'Q:\Research\fly_image\raw';
% target = 'Q:\Research\fly_image';
% splitp = 0.2;

rng(1) % for reproducability

% Get the folders in the raw data, which designate data labels
flisr = dir(root);
flisr = flisr(3:end);
n_group = size(flisr,1);

% Make train & test directories
traindir = fullfile(target, 'train');
testdir = fullfile(target, 'test');
mkdir(traindir)
mkdir(testdir)

% Make folders for each group for train & test data
for n = 1:n_group
    rawG = fullfile(root, flisr(n).name);
    trainG = fullfile(traindir, flisr(n).name);
    testG = fullfile(testdir, flisr(n).name);
    mkdir(trainG)
    mkdir(testG)
    
    imgG = dir(rawG);
    dI = cat(1, imgG.isdir);
    imgG = imgG(~dI);
    
    n_img = size(imgG,1);
    randI = randperm(n_img);
    n_test = floor(splitp*n_img);
    n_train = n_img - n_test;
    testI = randI(1:n_test);
    trainI = randI(n_test+1:end);
    
    for f = 1:n_test
        idx = testI(f);
        rawpath = fullfile(rawG, imgG(idx).name);
        newpath = fullfile(testG, imgG(idx).name);
        copyfile(rawpath, newpath);
    end
    
    for f = 1:n_train
        idx = trainI(f);
        rawpath = fullfile(rawG, imgG(idx).name);
        newpath = fullfile(trainG, imgG(idx).name);
        copyfile(rawpath, newpath);
    end
    
end
disp('Done')

end