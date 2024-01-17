function [data, mask] = abdomen_tracker(vid, save_path, mat_var_name,  mask_mode, npts, playback, vidpath)
%% abdomen_tracker: tracks insect abdomen yaw angle
%
% Tracks tip of abdomen, calculates the angle with
% respect to a specififed center point.
%
% Sign convention for angle outputs: [CW = + , CCW = -]
%
%   INPUT:
%       vid             : raw grayscale video data in matrix form
%                         or a path to a video file (.mp4 or .avi)
%       save_path       : path to save calculatd body angles.
%                         if empty, default is 'tracked_body' folder in 'vid' directory
%       mat_var_name    : name of the matlab variable containing the video data
%                         only required if inputing vid_path to .mat file
%       mask_mode 	    : [a b]. 'a' is a boolean that sets whether initial mask position is 
%                         automatically computed (true) or not (false). 'b' is a boolean that sets 
%                         whether initial mask position is editable by user after being set by 
%                         x (true) or left using the initial position (false). 
%                         Ex: [1 1] (set automatically & and let user edit after)
%                         Can also specifiy a mask structure saved from previous abdomen tracker.
%       npts            : # of points for tracker (if nan then track all points)
%       playback        : playback rate (show a frame in increments of "playback")
%                         if false, then don't show anything (default = 1)
%       vidpath         : path of output video. No ouput if false
%
%   OUTPUT:
%       data (structure)
%           angle       :   abdomen angles in body reference frame [°]
%           angle_glob	:   abdomen angles in global reference frame [°]
%           clust     	:   point cluster labels
%           points     	:   point locations
%           tip     	:   tracked tip points
%
%       mask: mask structure
%

% Parse inputs
if nargin < 7
    vidpath = [];
    if nargin < 6
        playback = 1;
        if nargin < 5
            npts = 50;
            if nargin < 4
                mask_mode = [1, 1];
                if nargin < 3
                    mat_var_name = 'vid';
                    if nargin < 2
                        save_path = [];
                    end
                end
            end
        end
    end
end

if isempty(mat_var_name)
    mat_var_name = 'vid';
end

% Parse video input
[viddata] = get_vid_data(vid, save_path, mat_var_name, 'tracked_abdomen');

% Set output video path
if ~isempty(vidpath) % export video to default path
    if ~vidpath % don't export
        export = false;
    else
        export = true;
    end
else
    export = true;
    vidpath = viddata.save_vid_path;
end

if ~rem(playback,1)==0
    warning('Warning: "playback" roundes to nearest integer')
    playback = round(playback);
end

if length(mask_mode) > 1
    auto_mask = mask_mode(2); % check manually based on input
else
    auto_mask = false; % check manually by default
end

dim = size(viddata.first_frame); % get size of video

% Set mask
warning('off', 'MATLAB:declareGlobalBeforeUse')
if isstruct(mask_mode(1)) % mask given
    mask = mask_mode;
elseif mask_mode(1) == 1 % set mask automatically by finding neck joint
 	%warning('auto abdomen mask not yet implemented, set mask_mod(1) = 0')
    
    rot_frame = imbinarize(mean(viddata.first_frame,3));
    [ry,rx] = find(rot_frame);
    pivot = [median(rx) , 1.2*median(ry)];
    R = max(sum(rot_frame,1));
    
    global mask
 	mask_frame = viddata.first_frame; % get 1st frame to set mask
    mask = make_mask(pivot, 180, [0.08*R 0.3*R], [50 50], mask_frame, [1 0 1]);
    
    if auto_mask % check manually
        mask.fig.main.Units = 'inches';
        mask.fig.main.Position(3:4) = [8, 8];
        movegui('center')
        uiwait(mask.fig.main)
    else
        close(mask.fig.main)
    end
elseif mask_mode(1) == 0 % set mask manually, start at center
    pivot = [dim(2), dim(1)] ./ 2; % default is center of frame
    R = 0.2*dim(1);
    
    global mask
	mask_frame = viddata.first_frame; % get 1st frame to set mask
    mask = make_mask(pivot, 180, [0.8*R 2.2*R], [40 40], mask_frame, [1 0 1]);
    mask.fig.main.Units = 'inches';
    mask.fig.main.Position(3:4) = [8, 8];
    movegui('center')
    uiwait(mask.fig.main) % wait to set
end
disp('Mask set')

% Track abdomen in all frames using tip tracker
tic
disp('Tracking...')
norm = 2;
n_clust = 1;
dthresh = 8;
rmv_out = true;
abdomen.angle = nan(viddata.n_frame,1);
abdomen.angle_glob = nan(viddata.n_frame,1);
abdomen.clust = cell(viddata.n_frame,1);
abdomen.points = cell(viddata.n_frame,1);
abdomen.tip = nan(viddata.n_frame,2);
pivot = mask.move_points.rot;
SE = strel("disk",5);
for n = 1:viddata.n_frame
    % Get frame
    if viddata.matflag % from matrix
        raw_frame = viddata.vid(:, :, n);
    else % from video reader 
        raw_frame = im2gray(read(viddata.Vread, n));
    end

    % Invert if backinglighting was detected
    if viddata.invertflag
        frame = imcomplement(raw_frame);
    else
        frame = raw_frame;
    end
    
    track_frame = imerode(frame, SE);

    [angle, ~, pts, k, ~, c_mean_all] = tracktip(track_frame, mask.area_points, ...
        pivot, norm, npts, 'clust', n_clust, dthresh, rmv_out);

    abdomen.angle_glob(n) = angle - 90;
    
    mask.global = rad2deg(wrapTo2Pi(deg2rad(mask.global)));
    
    abdomen.angle(n) = abdomen.angle_glob(n) + mask.global - 180;
    abdomen.points{n} = pts;
    abdomen.clust{n} = k;
    abdomen.tip(n,:) = c_mean_all;
end
abdomen.length = sqrt(sum((abdomen.tip - mask.move_points.rot).^2, 2));
toc

% Show tracking if set
if playback
    if export
        VID = VideoWriter(vidpath,'MPEG-4');
        VID.FrameRate = 60;
        VID.Quality = 100;
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
        ax(2) = subplot(3,4,9:12); hold on ; cla ; xlim([0 viddata.n_frame])
            xlabel('Frame')
            ylabel('Angle (°)')
            h.hAngle = animatedline(ax(2), 'Color', 'm', 'LineWidth', 1);
            % ylim(5*[-1 1])

    set(ax, 'Color', fColor, 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
        'YColor', aColor, 'XColor',aColor)
    set(fig, 'Visible', 'on')
    
    for n = 1:viddata.n_frame
        % Display frame count every every 100 frames
        if (n==1) || ~mod(n,100) || (n==viddata.n_frame)
            fprintf('%i\n',n)
        end

        % Get frame
        if viddata.matflag % from matrix
            raw_frame = viddata.vid(:, :, n);
        else % from video reader 
            raw_frame = im2gray(read(viddata.Vread, n));
        end

        % Display image
        if (n==1) || (~mod(n,playback)) || (n==viddata.n_frame) % at playback rate
            % Show images with tracking annotation
            ax(1) = subplot(3,4,1:8); cla % frame & tracking
                imshow(raw_frame) ; hold on ; %title(angle)
                patch(mask.points(:,1), mask.points(:,2), ...
                    mask.color, 'FaceAlpha', 0.2, 'EdgeColor', mask.color, 'LineWidth', 0.5);

                pts_left = abdomen.points{n}(abdomen.clust{n}==1,:);
                pts_right = abdomen.points{n}(abdomen.clust{n}==2,:);
                if isempty(pts_left) ||isempty(pts_right)
                    plot(abdomen.points{n}(:,1), abdomen.points{n}(:,2), ...
                        '.', 'MarkerSize', 5, 'Color', 0.5*mask.color)
                else
                    plot(pts_left(:,1), pts_left(:,2), 'r.', 'MarkerSize', 5)
                    plot(pts_right(:,1), pts_right(:,2), 'b.', 'MarkerSize', 5)
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
end
if export
    close(VID)
end

% Save abdomen data & mask
data = abdomen;
save(viddata.save_data_path, 'data', 'mask', '-v7.3')

end