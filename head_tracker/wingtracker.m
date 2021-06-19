function [wing, mask] = wingtracker(vid, mask_mode, playback)
% wingtracker: tracks insect head movments with two visible antenna
%
% Tracks antenna tips on fly head calculates the angle with
% respect to a specififed center point.
%
% Sign convention for angle outputs: [CW = + , CCW = -]
%
%   INPUT:
%       vid         :   video matrix
%       mask_mode 	:   [x y]. x is a boolean that sets whether initial mask position is 
%                       automatically computed (true) or not (false). y is a boolean that sets 
%                       whether initial mask position is editable by user after being set by 
%                       x (true) or left using the initial position (false). 
%                       Ex: [1 1] (set automatically & and let user edit after)
%                       Can also specifiy a mask structure saved from previous head tracker.
%       npts        :   # of points for tracker
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%
%   OUTPUT:
%       head (structure)
%           angle       :   head angles in body reference frame [°]
%           angle_glob	:   head angles in global reference frame [°]
%           tip     	:   tracked tip points
%
%       mask: mask structure
%

if nargin < 4
    playback = false;
    if nargin < 3
        npts = 10;
        if nargin < 2
            mask_mode = [1 1];
        end
    end
end

warning('off', 'MATLAB:declareGlobalBeforeUse')
if ~rem(playback,1)==0
    warning('Warning: "playback" roundes to nearest integer')
    playback = round(playback);
end

if length(mask_mode) > 1
    auto_mask = mask_mode(2); % check manually baed on input
else
    auto_mask = true; % don't check manually
end

vid = squeeze(vid); % remove singleton dimensions
dim = size(vid); % get size of video

% Set mask
if isstruct(mask_mode(1)) % mask given
    wing_mask = mask_mode;
elseif mask_mode(1) == 1 % set mask automatically by finding neck joint
 	disp('Finding neck ...')
    [VID,cent] = get_cut_vid(vid, 0.2, [], [0.1 0.1], 0.3); % get head vid
    [pivot,R,~,body_yaw] = get_neck(VID.out, VID.bw, cent); % find neck joint, head radius, & body angle
    
 	mask_frame = vid(:,:,1); % get 1st frame to set mask
    imshow(mask_frame)
    wing_mask = fly_mask(mask_frame, pivot, median(body_yaw));
    
    if auto_mask % check manually
        uiwait(gcf)
    else
        close(gcf)
    end
elseif mask_mode(1) == 0 % set mask manually, start at center
    pivot = [dim(2), dim(1)] ./ 2; % default is center of frame
    R = 0.1*dim(1);
    
	mask_frame = vid(:,:,1); % get 1st frame to set mask
    imshow(mask_frame)
    wing_mask = fly_mask(mask_frame, pivot, 0);
    uiwait(gcf) % wait to set
end
disp('Mask set')
% temp = draw(wing_mask);
% wing_mask = temp;
% close
% get_mask(wing_mask);
global mask

% Track in all frames using edge tracker
tic
disp('Tracking...')
norm = 2;
wing.angle = nan(dim(3),2);
wing.angle_glob = nan(dim(3),2);
wing.tip = nan(dim(3),4);
pivot = mask.lwing.move_points.rot;
for n = 1:dim(3)   
    [angle,tip] = trackedge(vid(:,:,n), mask.lwing.area_points, ...
        pivot, norm, false);
    wing.angle_glob(n) = angle - 270;
    wing.angle(n) = wing.angle_glob(n) - mask.global;
    wing.tip(n,1:2) = tip;
end
toc

% Show tracking if set
if playback
    % Create display
    fig(1) = figure (101); clf
    fColor = 'k'; % figure and main axis color
    aColor = 'w'; % axis line color
    set(fig, 'Color', fColor, 'Units', 'inches', 'Name', 'CrazyFly')
    fig.Position(3:4) = [9 7];
    movegui(fig, 'center')
    figure (101)
        % Raw image with tracking window
        ax(1) = subplot(3,4,1:8); hold on ; cla ; axis image

        % Tracking angle window
        ax(2) = subplot(3,4,9:12); hold on ; cla ; xlim([0 dim(3)])
            xlabel('Frame')
            ylabel('Angle (°)')
            h.hAngle = animatedline(ax(2), 'Color', 'c', 'LineWidth', 1);

    set(ax, 'Color', fColor, 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
        'YColor', aColor, 'XColor',aColor)
    set(fig, 'Visible', 'on')
    
    for n = 1:dim(3)
        % Display frame count every every 100 frames
        if (n==1) || ~mod(n,100) || (n==dim(3))
            fprintf('%i\n',n)
        end

        % Display image
        if (n==1) || (~mod(n,playback)) || (n==dim(3)) % at playback rate
            % Show images with tracking annotation
            ax(1) = subplot(3,4,1:8); cla % frame & tracking
                imshow(imadjust(vid(:,:,n))) ; hold on ; title(angle)
                patch(mask.points(:,1), mask.points(:,2), ...
                    mask.color, 'FaceAlpha', 0.2, 'EdgeColor', mask.color, 'LineWidth', 0.5);

                plot([pivot(1) mask.move_points.axis(1)], [pivot(2) mask.move_points.axis(2)], ...
                    '--', 'Color', 0.5*mask.color, 'LineWidth', 0.5)
                
                plot([pivot(1) wing.tip(n,1)], [pivot(2) wing.tip(n,2)], ...
                    'Color', mask.color, 'LineWidth', 2)

                plot(pivot(1), pivot(2), '.', 'Color', mask.color, 'MarkerSize', 10)
        end

        % Display angle
        ax(2) = subplot(3,4,9:12); % angles
            addpoints(h.hAngle, n, wing.angle(n))

        pause(0.0005) % give time for images to display
    end
end

end