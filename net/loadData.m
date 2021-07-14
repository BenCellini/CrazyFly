function [train,test,val] = loadData(dataDir,train_folder,test_folder)
%LOADDATA load training and testing data sets
%   

rng(1) % For reproducibility

Groups = {'Up', 'Down'};

fprintf('Loading Train Filenames and Label Data...'); t = tic;
train_all = imageDatastore(fullfile(dataDir,train_folder),...
    'IncludeSubfolders',true,'LabelSource','foldernames');
train_all.Labels = reordercats(train_all.Labels,Groups);

% Split with validation set
[train, val] = splitEachLabel(train_all,.9);
fprintf('Done in %.02f seconds\n', toc(t));

fprintf('Loading Test Filenames and Label Data...'); t = tic;
test = imageDatastore(fullfile(dataDir,test_folder),...
    'IncludeSubfolders',true,'LabelSource','foldernames');
test.Labels = reordercats(test.Labels,Groups);
fprintf('Done in %.02f seconds\n', toc(t));
end

