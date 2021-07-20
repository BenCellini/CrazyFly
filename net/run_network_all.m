%% Normal Dataset
clear
close all
clc

dataDir = 'Q:\Research\fly_image_frye\';
train_folder = 'train';
test_folder  = 'test';
log = fullfile(dataDir, 'log');

% Load the data
[train,test,val] = loadData(dataDir, train_folder, test_folder);

% Set CNN parameters
numEpochs = 12;
batchSize = 50;
learnRate = 1e-6;

% batchSize = fliplr(50:25:500);
% learnRate = fliplr([1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4]);

% Train network
ALL = cell(length(batchSize),length(learnRate));
for b = 1:length(batchSize)
    for r = 1:length(learnRate)
        [status, message, messageid] = rmdir('modelCheckpoints', 's');
        ALL{b,r} = trainNormal(train, test, val, numEpochs, batchSize(b), learnRate(r), log);
        close all
    end
end

%% 
clear ; close all ; clc
dataDir = 'Q:\Research\fly_image\';
train_folder = 'train';
test_folder  = 'test';
log = fullfile(dataDir, 'log');

% Load the data
[train,test,val] = loadData(dataDir, train_folder, test_folder);

load('Q:\Research\fly_image\log\Normal_numEpoch_12_batchSize_50_learnRate_1e-06.mat', 'network')

% Y_test = classify(network,test);
% test_acc = mean(Y_test==test.Labels);
I_test = imread("Q:\Research\fly_image\test\Up\Experiment_SS_vel_250_fly_1_trial_6_amp_60_freq_0.7_frame_628.jpg");
Y_test = classify(network,I_test);

%%
close all
clc

scale_nn = 1;
NN = load('Q:\Research\fly_image_frye\log\Normal_numEpoch_12_batchSize_50_learnRate_1e-06.mat', ...
    'network', 'ALL');
nn_sz = NN.network.Layers(1).InputSize(1:2);

yx = ceil(scale_nn*nn_sz);

vid = squeeze(vidData);
dim = size(vid);
test_up = uint8(zeros(yx));
bAngles = nan(dim(3),1);
tic
for n = 1:dim(3)
    frame = vid(:,:,n);
    [heading,fly_frame] = getflyroi_ud(frame, yx);
    input_frame = imresize(fly_frame, nn_sz);
    Y = classify(NN.network, input_frame);
    switch Y
        case 'Up'
            out_frame = fly_frame;
        case 'Down'
            out_frame = rot90(fly_frame,2);
            heading = heading + 180;
    end
    bAngles(n) = heading;
    test_up(:,:,n) = imresize(out_frame, yx);
end
toc

