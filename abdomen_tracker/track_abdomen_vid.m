function [abdomen, mask] = track_abdomen_vid(vid, set_mask, npts, playback, vidpath)
% track_abdomen: tracks insect abdomen movements
%
% Tracks tip of abdomen, calculates the angle with
% respect to a specififed center point.
%
% Sign convention for angle outputs: [CW = + , CCW = -]
%
%   INPUT:
%       vid         :   video matrix
%       set_mask  	:   predefined mask or 2x1: [1 1] = auto place mask and let user adjust
%       npts        :   # of points for tracker (if nan then track all points)
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%
%   OUTPUT:
%       abdomen (structure)
%           angle       :   abdomen angles in body reference frame [°]
%           angle_glob	:   abdomen angles in global reference frame [°]
%           clust     	:   point cluster labels
%           points     	:   point locations
%           tip     	:   tracked tip points
%
%       mask: mask structure
%

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

if ~isempty(vidpath)
    export = true;
else
    export = false;
end

vid = squeeze(vid); % remove singleton dimensions
dim = size(vid); % get size of video

% Set mask
if isstruct(set_mask(1)) % mask given
    mask = set_mask;
elseif set_mask(1) == 1 % set mask automatically by finding neck joint
    %error('this ''set_mask'' method doesnt exist yet')
 	%disp('Finding neck ...')
    
    rot_frame = imbinarize(mean(vid,3));
    [ry,rx] = find(rot_frame);
    pivot = [median(rx) , 1.2*median(ry)];
    R = max(sum(rot_frame,1));
    
    global mask
 	mask_frame = vid(:,:,1); % get 1st frame to set mask
    mask = make_mask(pivot, 180, [0.08*R 0.3*R], [50 50], mask_frame, [1 0 1]);
    
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
    mask = make_mask(pivot, 180, [0.8*R 2.2*R], [40 40], mask_frame, [1 0 1]);
    uiwait(mask.fig.main) % wait to set
end
disp('Mask set')
% close

% Track abdomen in all frames using tip tracker
tic
disp('Tracking...')
norm = 2;
n_clust = 1;
dthresh = 8;
rmv_out = true;
abdomen.angle = nan(dim(3),1);
abdomen.angle_glob = nan(dim(3),1);
% abdomen.antenna = nan(dim(3),n_clust);
abdomen.clust = cell(dim(3),1);
abdomen.points = cell(dim(3),1);
abdomen.tip = nan(dim(3),2);
pivot = mask.move_points.rot;
SE = strel("disk",5);
for n = 1:dim(3)
    frame = vid(:,:,n);
    track_frame = imerode(frame,SE);
    %imshow(track_frame)
    [angle,m,pts,k] = tracktip(track_frame, mask.area_points, ...
        pivot, norm, npts, 'clust', n_clust, dthresh, rmv_out);
    abdomen.angle_glob(n) = angle - 90 + 360;
    abdomen.angle(n) = abdomen.angle_glob(n) - mask.global - 180;
    abdomen.points{n} = pts;
    abdomen.clust{n} = k;

    abdomen.tip(n,:) = mask.move_points.rot +  ...
        -mask.radius.outer*[sind(abdomen.angle_glob(n)) , -cosd(abdomen.angle_glob(n))];
end
toc

% Show tracking if set
if playback
    if export
        VID = VideoWriter(vidpath,'MPEG-4');
        VID.FrameRate = 100;
        open(VID)
    end
    
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
        ax(2) = subplot(3,4,9:12); hold on ; cla ; xlim([0 dim(3)])
            xlabel('Frame')
            ylabel('Angle (°)')
            h.hAngle = animatedline(ax(2), 'Color', 'm', 'LineWidth', 1);
            % ylim(5*[-1 1])

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
                imshow(vid(:,:,n)) ; hold on ; title(angle)
                patch(mask.points(:,1), mask.points(:,2), ...
                    mask.color, 'FaceAlpha', 0.2, 'EdgeColor', mask.color, 'LineWidth', 0.5);

                pts_left = abdomen.points{n}(abdomen.clust{n}==1,:);
                pts_right = abdomen.points{n}(abdomen.clust{n}==2,:);
                if isempty(pts_left) ||isempty(pts_right)
                    plot(abdomen.points{n}(:,1),abdomen.points{n}(:,2), 'm.', 'MarkerSize', 5)
                else
                    plot(pts_left(:,1),pts_left(:,2), 'r.', 'MarkerSize', 5)
                    plot(pts_right(:,1),pts_right(:,2), 'b.', 'MarkerSize', 5)
                end

                plot([pivot(1) mask.move_points.axis(1)], [pivot(2) mask.move_points.axis(2)], ...
                    '--', 'Color', 0.5*mask.color, 'LineWidth', 0.5)
                
                plot([pivot(1) abdomen.tip(n,1)], [pivot(2) abdomen.tip(n,2)], ...
                    'Color', mask.color, 'LineWidth', 2)

                plot(pivot(1), pivot(2), '.', 'Color', mask.color, 'MarkerSize', 10)
        end

        % Display angle
        ax(2) = subplot(3,4,9:12); % angles
            addpoints(h.hAngle, n, abdomen.angle(n))

        drawnow % give time for images to display
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