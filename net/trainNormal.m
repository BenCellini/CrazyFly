function [ALL] = trainNormal(train, test, val, numEpochs, batchSize, learnRate, logdir)
%% trainNormal: Train a CNN on the unaugmented dataset
%   Trains a CNN
%
%   INPUT:
%       train       : training data
%       test        : testing data
%       numEpochs   : # of epochs
%       batchSize   : batch size
%       learnRate   : learning rate
%
%   OUTPUT:
%       ALL        	: structure containing netwrok properties
%

rng('default');

% Make the checkpoint directory
checkpointDir = 'modelCheckpoints';
if ~exist(checkpointDir,'dir'); mkdir(checkpointDir); end

nTraining = length(train.Labels);
     
% Define the Network Structure, To add more layers, copy and paste the
% lines such as the example at the bottom of the code
% CONV -> ReLU -> POOL -> FC -> DROPOUT -> FC -> SOFTMAX 
layers = [
    imageInputLayer([222 133 1]); % Input to the network is a 256x256x1 sized image 
    convolution2dLayer(5,20,'Padding',[2 2],'Stride', [2,2]);  % convolution layer with 20, 5x5 filters
    reluLayer(); % ReLU layer
    maxPooling2dLayer(2,'Stride',2); % Max pooling layer
    fullyConnectedLayer(25); % Fullly connected layer with 50 activations
    dropoutLayer(0.25); % Dropout layer
    fullyConnectedLayer(2); % Fully connected with 17 layers
    softmaxLayer(); % Softmax normalization layer
    classificationLayer(); % Classification layer
    ];

% Set the training options
options = trainingOptions('sgdm','InitialLearnRate', learnRate,...% learning rate
    'CheckpointPath', checkpointDir,...
    'MiniBatchSize', batchSize, ...
    'MaxEpochs',numEpochs,...
    'ExecutionEnvironment','auto', ...
    'OutputFcn',@plotTrainingAccuracy);

% Train the network, info contains information about the training accuracy and loss
t = tic;
[network,info] = trainNetwork(train,layers,options);
trainTime = toc(t);
fprintf('Trained in in %.02f seconds\n', trainTime);

% Test on the validation data
Y_train = classify(network,val);
val_acc = mean(Y_train==val.Labels);
fprintf('Training Accuracy: %f \n', val_acc);

% Test on the Testing data
Y_test = classify(network,test);
test_acc = mean(Y_test==test.Labels);

% Make structure to store network information
ALL.network     = network;
ALL.info        = info;
ALL.layers      = layers;
ALL.options     = options;
ALL.trainTime   = trainTime;
ALL.Y.train     = Y_train;
ALL.Y.test      = Y_test;
ALL.val_acc     = val_acc;
ALL.test_acc   	= test_acc;

% Plot
fig = figure (101); clf
plotTrainingAccuracy_All(info,numEpochs);

% Save figure and data
fname = ['Normal_numEpoch_' num2str(numEpochs) '_batchSize_' num2str(batchSize) ... % filename
            '_learnRate_' num2str(learnRate)];
     
if ~isempty(logdir)
    data_file = fullfile(logdir,[fname '.mat']);
    fig_file  = fullfile('figure',[fname '.fig']);
    save(data_file, 'ALL', 'network', 'info', 'numEpochs', 'batchSize', 'learnRate', 'nTraining', 'trainTime', ...
                'Y_train', 'Y_test', 'val_acc', 'test_acc', 'layers', 'options', 'fname');
    savefig(fig, fig_file)
end

end

