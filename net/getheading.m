function [body_ang, stable_vid] = getheading(vid, scale, showplot, NN)
%% get_heading: uses a trained neural network to find the heading of a fly in an video/image
%
%   INPUT:
%       vid         : video or image
%       scale       : scale the image ROI by this factor if the image around the fly is too large/small
%       showplot   	: show heading tracking (boolean)
%       NN          : trained neural network, setting this input should save computation time
%
%   OUTPUT:
%       body_ang    : body angle [deg]
%       stable_vid  : stabilized video
%

if nargin < 4
    %----------------------------------------------------------------------------------------------
    NN = load('Q:\Research\fly_image_frye\log\Normal_numEpoch_12_batchSize_50_learnRate_1e-06.mat', ...
        'network'); % EDIT THIS TO LOCATION ON YOUR PC
    %----------------------------------------------------------------------------------------------
    if nargin < 3
        showplot = false;
        if nargin < 2
            scale = 1;
        end
    end
end

tic

% Get input image size for NN
nn_sz = NN.network.Layers(1).InputSize(1:2);
yx = ceil(scale*nn_sz);

% Format input video
vid = squeeze(vid);
dim = size(vid);
if length(dim) < 3
   dim(3) = 1; 
end

% Set up variables to store outputs
stable_vid = uint8(zeros(yx));
body_ang = nan(dim(3),1);
if showplot
    fig = figure (111);
    ax(1) = subplot(2,2,1); cla; hold on; title('raw')
    ax(2) = subplot(2,2,2); cla; hold on; title('stable')
    ax(2) = subplot(2,2,3:4); cla; hold on; title('body angle (deg)')
        h.body = animatedline('Color', 'r', 'LineWidth', 1.5);
end

% Find heading in each frame
for n = 1:dim(3)
    frame = vid(:,:,n);
    [heading, fly_frame, imgstats] = getflyroi_ud(frame, yx); % get the rotated fly frame
    input_frame = imresize(fly_frame, nn_sz); % resize to fit NN
    Y = classify(NN.network, input_frame); % classify using NN
    
    % Flip heading by 180 deg if fly is upside down
    switch Y
        case 'Up'
            out_frame = fly_frame;
        case 'Down'
            out_frame = rot90(fly_frame,2);
            heading = heading + 180;
    end
    body_ang(n) = heading;
    stable_vid(:,:,n) = imresize(out_frame, yx);
    
    % Show figure
    if showplot
        cent = imgstats.Centroid;
        R = imgstats.MajorAxisLength / 4;
        head_tip = cent + R*[cosd(body_ang(n)) , -sind(body_ang(n))];
        subplot(2,2,1)
            imshow(vid(:,:,n))
            plot([cent(1) head_tip(1)], [cent(2) head_tip(2)], 'r', 'LineWidth', 1)
            plot(cent(1), cent(2), 'r.', 'MarkerSize', 15)
        subplot(2,2,2)
            imshow(stable_vid(:,:,n))
        subplot(2,2,3:4)
            addpoints(h.body, n, body_ang(n))
            drawnow
    end    
end
toc

end