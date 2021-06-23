function [bodypart, mask, offsetAngle] = track_any(vid, set_mask, npts, offsetAngle, playback, vidpath)
% track_any: tracks insect body part
%
% Tracks tip of bodypart, calculates the angle with
% respect to a specififed center point.
%
% Sign convention for angle outputs: [CW = + , CCW = -]
%
%   INPUT:
%       vid         :   video matrix
%       set_mask  	:   predefined mask or 2x1: [1 1] = auto place mask and let user adjust
%       npts        :   # of points for tracker (if nan then track all points)
%       offsetAngle :   offset angle from tracked tip to true angle
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%       export      : file path to export to .mp4
%
%   OUTPUT:
%       bodypart (structure)
%           angle       :   angles in body reference frame [°]
%           angle_glob	:   angles in global reference frame [°]
%           clust     	:   point cluster labels
%           points     	:   point locations
%           tip     	:   tracked tip points
%
%       mask: mask structure
%

if nargin < 6
    vidpath = [];
    if nargin < 5
        playback = true;
        if nargin < 4
            offsetAngle = 0;
        end
    end
end

if ~isempty(vidpath)
    export = true;
else
    export = false;
end

warning('off', 'MATLAB:declareGlobalBeforeUse')
if ~rem(playback,1)==0
    warning('Warning: "playback" roundes to nearest integer')
    playback = round(playback);
end

if length(set_mask) > 1 % check manually
    auto_mask = set_mask(2);
else
    auto_mask = true; 
end

vid = squeeze(vid); % remove singleton dimensions
dim = size(vid); % get size of video

% Set mask
if isstruct(set_mask(1)) % mask given
    mask = set_mask;
elseif set_mask(1) == 1 % set mask automatically by finding neck joint    
    rot_frame = imbinarize(mean(vid,3));
    [ry,rx] = find(rot_frame);
    pivot = [median(rx) , 1.2*median(ry)];
    R = max(sum(rot_frame,1));
    
    global mask
 	mask_frame = vid(:,:,1); % get 1st frame to set mask
    mask = make_mask(pivot, 0, [0.08*R 0.3*R], [50 50], mask_frame, [0 1 0.5]);
    
    if auto_mask % check manually
        uiwait(mask.fig.main)
    else
        close(mask.fig.main)
    end
elseif set_mask(1) == 0 % set mask manually, start at center
    pivot = [dim(2), dim(1)] ./ 2; % default is center of frame
    R = 0.2*dim(1);
    
    global mask
	mask_frame = vid(:,:,1); % get 1st frame to set mask
    mask = make_mask(pivot, 0, [0.8*R 2.2*R], [40 40], mask_frame, [0 1 0.5]);
    uiwait(mask.fig.main) % wait to set
end
disp('Mask set')
% close

% Track bodypart in all frames using tip tracker
tic
disp('Tracking...')
norm = 2;
n_clust = 1;
dthresh = 8;
rmv_out = true;
bodypart.angle = nan(dim(3),1);
bodypart.angle_glob = nan(dim(3),1);
% abdomen.antenna = nan(dim(3),n_clust);
bodypart.clust = cell(dim(3),1);
bodypart.points = cell(dim(3),1);
bodypart.tip = nan(dim(3),2);
pivot = mask.move_points.rot;
% SE = strel("disk",5);
for n = 1:dim(3)
    frame = vid(:,:,n);
    track_frame = frame;
    %track_frame = imerode(frame,SE);
    %imshow(track_frame)
    [angle,m,pts,k] = tracktip(track_frame, mask.area_points, ...
        pivot, norm, npts, 'clust', n_clust, dthresh, rmv_out);
    bodypart.points{n} = pts;
    bodypart.clust{n} = k;
    if isempty(offsetAngle) && (n == 1)
       offsetAngle = angle - 270 - mask.global - mask.init.angle; 
    end
    bodypart.angle_glob(n) = angle - 270 - offsetAngle;
    bodypart.angle(n) = bodypart.angle_glob(n) - mask.global;
    
    bodypart.tip(n,:) = mask.move_points.rot +  ...
        mask.radius.outer*[sind(bodypart.angle_glob(n)) , -cosd(bodypart.angle_glob(n))];
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

        % Head angle window
        ax(2) = subplot(3,4,9:12); hold on ; cla ; xlim([-0.02*dim(3) dim(3)])
            xlabel('Frame')
            ylabel('Angle (°)')
            h.hAngle = animatedline(ax(2), 'Color', mask.color, 'LineWidth', 1);
            % ylim(5*[-1 1])

    set(ax, 'Color', fColor, 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
        'YColor', aColor, 'XColor',aColor)
    set(fig, 'Visible', 'on')
    
    if export
        VID = VideoWriter(vidpath,'MPEG-4');
        VID.FrameRate = 100;
        open(VID)
    end
    
    for n = 1:dim(3)
        % Display frame count every every 100 frames
        if (n==1) || ~mod(n,100) || (n==dim(3))
            fprintf('%i\n',n)
        end

        % Display image
        if (n==1) || (~mod(n,playback)) || (n==dim(3)) % at playback rate
            % Show images with tracking annotation
            ax(1) = subplot(3,4,1:8); cla % frame & tracking
                imshow(vid(:,:,n)) ; hold on ; title(angle)
                patch(mask.points(:,1), mask.points(:,2), ...
                    mask.color, 'FaceAlpha', 0.2, 'EdgeColor', mask.color, 'LineWidth', 0.5);

                pts_left = bodypart.points{n}(bodypart.clust{n}==1,:);
                pts_right = bodypart.points{n}(bodypart.clust{n}==2,:);
                if isempty(pts_left) ||isempty(pts_right)
                    plot(bodypart.points{n}(:,1),bodypart.points{n}(:,2), 'b.', 'MarkerSize', 5)
                else
                    plot(pts_left(:,1),pts_left(:,2), 'r.', 'MarkerSize', 5)
                    plot(pts_right(:,1),pts_right(:,2), 'b.', 'MarkerSize', 5)
                end

                plot([pivot(1) mask.move_points.axis(1)], [pivot(2) mask.move_points.axis(2)], ...
                    '--', 'Color', 0.5*mask.color, 'LineWidth', 0.5)
                
                plot([pivot(1) bodypart.tip(n,1)], [pivot(2) bodypart.tip(n,2)], ...
                    'Color', mask.color, 'LineWidth', 2)

                plot(pivot(1), pivot(2), '.', 'Color', mask.color, 'MarkerSize', 10)
        end

        % Display angle
        ax(2) = subplot(3,4,9:12); % angles
            addpoints(h.hAngle, n, bodypart.angle(n))

        pause(0.0005) % give time for images to display
        
        if export
          	fig_frame = getframe(fig);
         	writeVideo(VID, fig_frame); 
        end
    end
    if export
       close(VID) 
    end
end

end