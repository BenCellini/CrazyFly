function [data, mask] = head_tracker(vid, save_path, mat_var_name, mask_mode, npts, neck_frames, playback, vidpath)
% track_head_vid: tracks head movments
%
% Tracks antenna tips on fly head calculates the angle with
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
%                         Can also specifiy a mask structure saved from previous head tracker.
%       npts            : # of points for tracker
%       neck_frames     : frames to use to find neck if mask_mode(1) = 1.
%                         if length(neck_frames) = 1, then choose that # of frames
%       playback        : playback rate (show a frame in increments of "playback")
%                         if false, then don't show anything (default = 1)
%       vidpath         : path of output video. No ouput if false
%
%   OUTPUT:
%       data (structure)
%           angle       :   head angles in body reference frame [deg]
%           angle_glob	:   head angles in global reference frame [deg]
%           antenna     :   antenna positions
%           clust     	:   point cluster labels
%           points     	:   point locations
%           tip     	:   tracked tip points
%
%       mask: mask structure
%

% Parse inputs
if nargin < 8
    vidpath = [];
    if nargin < 7
        playback = 1;
        if nargin < 6
            neck_frames = 20;
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
end

if isempty(mat_var_name)
    mat_var_name = 'vid';
end

% Parse video input
[viddata] = get_vid_data(vid, save_path, mat_var_name, 'tracked_head');

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

% Get frames to find neck
if viddata.matflag
    if length(neck_frames) == 1
        step = round(viddata.n_frame / neck_frames);
        neck_vid = viddata.vid(:, :, 1:step:viddata.n_frame);
    else
        neck_vid = viddata.vid(:, :, neck_frames);
    end
else
    if length(neck_frames) == 1
        step = round(viddata.n_frame / neck_frames);
        neck_vid = uint8(zeros([size(viddata.first_frame), neck_frames]));
        neck_frames = 1:step:viddata.n_frame;
        for f = 1:length(neck_frames)
            neck_vid(:, :, f) = im2gray(read(viddata.Vread, neck_frames(f)));
        end
    else
        neck_vid = uint8(zeros([size(viddata.first_frame), length(neck_frames)]));
        for f = neck_frames
            neck_vid(:, :, f) = im2gray(read(viddata.Vread, f));
        end
    end
end

% Set mask
warning('off', 'MATLAB:declareGlobalBeforeUse')
if isstruct(mask_mode(1)) % mask given directly
    mask = mask_mode;

elseif mask_mode(1) == 1 % set mask automatically by finding neck joint
 	disp('Finding neck ...')
    [VID,cent] = get_cut_vid(neck_vid, 1.0, [], [0.1 0.1], 0.3); % get head vid
    [pivot, R, ~, body_yaw] = get_neck(VID.out, VID.bw, cent); % find neck joint, head radius, & body angle
    
    global mask
 	mask_frame = viddata.first_frame; % get 1st frame to set mask
    mask = make_mask(pivot, median(body_yaw), [0.8*R 1.7*R], [40 40], mask_frame);
    mask.fig.main.Name = 'CrazyFly: Head-Tracker';
    
    if auto_mask % check manually
        mask.fig.main.Units = 'inches';
        mask.fig.main.Position(3:4) = [8, 8];
        text(0,-20,'Set the mask & click Done.', ...
            'Color', 'R', 'FontSize', 15, 'FontWeight', 'bold');
        movegui('center')
        uiwait(mask.fig.main)
    else
        close(mask.fig.main)
    end

elseif mask_mode(1) == 0 % set mask manually, start at center
    pivot = [dim(2), dim(1)] ./ 2; % default is center of frame
    R = 0.1*dim(1);
    
    global mask
	mask_frame = viddata.first_frame; % get 1st frame to set mask
    mask = make_mask(pivot, 0, [0.8*R 2.0*R], [40 40], mask_frame);
    mask.fig.main.Units = 'inches';
    mask.fig.main.Position(3:4) = [8, 8];
        text(0,-20,'Set the mask & click Done.', ...
            'Color', 'R', 'FontSize', 15, 'FontWeight', 'bold');
    movegui('center')
    uiwait(mask.fig.main) % wait to set
end
disp('Mask set')

% Track head in all frames using tip tracker
tic
disp('Tracking...')
norm = 2;
head.angle = nan(viddata.n_frame, 1);
head.angle_glob = nan(viddata.n_frame, 1);
head.antenna = nan(viddata.n_frame, 2);
head.clust = cell(viddata.n_frame, 1);
head.points = cell(viddata.n_frame, 1);
head.tip = nan(viddata.n_frame, 2);
pivot = mask.move_points.rot;
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

    % Run tip-tracker
    [angle, m, pts, k, ~, c_mean_all] = ...
    tracktip(frame, mask.area_points, pivot, norm, npts, 'clust', 2, 29);

    head.angle_glob(n) = angle - 270;
    head.angle(n) = head.angle_glob(n) - mask.global;
    head.antenna(n, :) = m' - 270;
    head.points{n} = pts;
    head.clust{n} = k;

    head.tip(n, :) = c_mean_all;

    %head.tip(n,:) = mask.move_points.rot +  ...
        %mask.radius.outer*[sind(head.angle_glob(n)) , -cosd(head.angle_glob(n))];
end
toc

% Show tracking if set
if playback
    if export
        VID = VideoWriter(vidpath,'MPEG-4');
        VID.Quality = 100;
        VID.FrameRate = 60;
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
            ylabel('Angle (\circ)', 'Interpreter', 'tex')
            h.hAngle = animatedline(ax(2), 'Color', 'c', 'LineWidth', 1);
            % ylim(5*[-1 1])

    set(ax, 'Color', fColor, 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
        'YColor', aColor, 'XColor',aColor)
    set(fig, 'Visible', 'on')
    
    for n = 1:viddata.n_frame
        % Display frame count every every 100 frames
        if (n==1) || ~mod(n,100) || (n == viddata.n_frame)
            fprintf('%i\n',n)
        end

        % Get frame
        if viddata.matflag % from matrix
            raw_frame = viddata.vid(:, :, n);
        else % from video reader 
            raw_frame = im2gray(read(viddata.Vread, n));
        end

        % Display image
        if (n==1) || (~mod(n,playback)) || (n == viddata.n_frame) % at playback rate
            % Show images with tracking annotation
            ax(1) = subplot(3,4,1:8); cla % frame & tracking
                imshow(raw_frame) ; hold on ; title(angle)
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

                plot([pivot(1) mask.move_points.axis(1)], [pivot(2) mask.move_points.axis(2)], ...
                    '--', 'Color', 0.5*mask.color, 'LineWidth', 0.5)
                
                plot([pivot(1) head.tip(n,1)], [pivot(2) head.tip(n,2)], ...
                    'Color', mask.color, 'LineWidth', 2)

                plot(pivot(1), pivot(2), '.', 'Color', mask.color, 'MarkerSize', 10)
        end

        % Display angle
        ax(2) = subplot(3,4,9:12); % angles
            addpoints(h.hAngle, n, head.angle(n))
            
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

% Save head data & mask
data = head;
save(viddata.save_data_path, 'data', 'mask', '-v7.3')

end