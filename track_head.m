function [head, mask] = track_head(vid, set_mask, npts, playback)
% track_head: tracks insect head movments
%
% Tracks antenna tips on fly head calculates the angle with
% respect to a specififed center point.
%
% Sign convention for angle outputs: [CW = + , CCW = -]
%
%   INPUT:
%       vid         :   video matrix
%       mask      	:   peredefined mask
%       npts        :   # of points for tracker
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%
%   OUTPUT:
%       head (structure)
%           angle       :   head angles [°]
%           mask        :   head mask
%           antenna   	:   antenna angles [°]
%           points     	:   tracked tip points
%           clust     	:   point cluster labels
%

warning('off', 'MATLAB:declareGlobalBeforeUse')
if ~rem(playback,1)==0
    warning('Warning: "playback" roundes to nearest integer')
    playback = round(playback);
end

if length(set_mask) > 1  % check manually
    auto_mask = set_mask(2);
else
    auto_mask = true; 
end

vid = squeeze(vid); % remove singleton dimensions
dim = size(vid); % get size of video

% Set mask
if isstruct(set_mask(1)) % mask given
    mask = set_mask;
    % no computations needed
elseif set_mask(1) == 1 % set mask automatically by finding neck joint
 	disp('Finding neck ...')
    %[body_yaw, ~, ~] = get_vid_props(vid, false); % bouding box
    [VID,cent] = get_cut_vid(vid, 0.2, [], [0.1 0.1], 0.3); % get head vid
    [pivot,R,~,body_yaw] = get_neck(VID.out, VID.bw, cent); % find neck joint % head radius
    
    global mask
 	mask_frame = vid(:,:,1); % get 1st frame to set mask
    %body_yaw = 0;
    mask = make_mask(pivot, median(body_yaw), [0.8*R 2.2*R], [40 40], mask_frame);
    
    if auto_mask % check manually
        uiwait(mask.fig.main)
    else
        close(mask.fig.main)
    end
elseif set_mask(1) == 2 % set mask manually, start at center
    pivot = [dim(2), dim(1)] ./ 2; % default is center of frame
    R = 0.1*dim(1);
    
    global mask
	mask_frame = vid(:,:,1); % get 1st frame to set mask
    mask = make_mask(pivot, 0, [0.8*R 2.2*R], [40 40], mask_frame);
    uiwait(mask.fig.main) % wait to set
end
disp('Mask set')

% Track head in all frames using tip tracker
tic
disp('Tracking...')
norm = 2;
head.angle = nan(dim(3),1);
head.angle_glob = nan(dim(3),1);
head.antenna = nan(dim(3),2);
head.clust = cell(dim(3),1);
head.points = cell(dim(3),1);
head.tip = nan(dim(3),2);
pivot = mask.move_points.rot;
for n = 1:dim(3)   
    [angle,m,pts_left,k] = tracktip(vid(:,:,n), mask.area_points, ...
        pivot, norm, npts, 'clust');
    head.angle_glob(n) = angle - 270;
    head.angle(n) = head.angle_glob(n) - mask.global;
    head.antenna(n,:) = m' - 270;
    head.points{n} = pts_left;
    head.clust{n} = k;

    head.tip(n,:) = mask.move_points.rot +  ...
        mask.radius.outer*[sind(head.angle(n)) , -cosd(head.angle(n))];
end
toc

% Show tracking if set
if playback
    % Create display
    fig(1) = figure (101); clf
    fColor = 'k'; % figure and main axis color
    aColor = 'w'; % axis line color
    set(fig,'Color',fColor,'Units','inches','Position',[2 2 9 7])
    figure (101)
        % Raw image with tracking window
        ax(1) = subplot(3,4,1:8); hold on ; cla ; axis image

        % Head angle window
        ax(2) = subplot(3,4,9:12); hold on ; cla ; xlim([0 dim(3)])
            xlabel('Frame')
            ylabel('Angle (°)')
            h.hAngle = animatedline(ax(2), 'Color', 'c', 'LineWidth', 1.5);
            % ylim(5*[-1 1])

    set(ax, 'Color', fColor, 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
        'YColor', aColor, 'XColor',aColor)
    set(fig, 'Visible', 'on')
    
    for n = 1:dim(3)
        % Display frame count every every 100 frames
        if (n==1) || ~mod(n,100) || (n==dim(3))
            fprintf('%i\n',n)
        end

        if playback || n==1 || n==dim(3)

            % Display image
            if (n==1) || (~mod(n,playback)) || (n==dim(3)) % at playback rate
                % Show images with tracking annotation
                ax(1) = subplot(3,4,1:8); cla % frame & tracking
                    imshow(vid(:,:,n)) ; hold on ; title(angle)
                    patch(mask.points(:,1), mask.points(:,2), ...
                        mask.color, 'FaceAlpha', 0.2, 'EdgeColor', mask.color, 'LineWidth', 0.5);
                    
                    pts_left = head.points{n}(head.clust{n}==1,:);
                    pts_right = head.points{n}(head.clust{n}==2,:);
                    if isempty(pts_left) ||isempty(pts_right)
                        plot(head.points{n}(:,1),head.points{n}(:,2), 'm.', 'MarkerSize', 5)
                    else
                        plot(pts_left(:,1),pts_left(:,2), 'r.', 'MarkerSize', 5)
                        plot(pts_right(:,1),pts_right(:,2), 'b.', 'MarkerSize', 5)
                    end
                    
                    plot([mask.move_points.rot(1) head.tip(n,1)], ...
                        [mask.move_points.rot(2) head.tip(n,2)], ...
                        'Color', mask.color, 'LineWidth', 2)

                    plot(mask.move_points.rot(1), mask.move_points.rot(2), '.', ...
                        'Color', mask.color, 'MarkerSize', 10)
            end

            % Display angle
            ax(2) = subplot(3,4,9:12); % angles
                addpoints(h.hAngle, n, head.angle(n))

            pause(0.0005) % give time for images to display
        end
    end
end

% Get stabilized head window
%[stable.vid, stable.cent] = stable_head(vid, head.angle, pivot);

end