function [angle, imgstats, initframe] = body_tracker(vid, save_path, mat_var_name, debug_heading, flip_vid, playback, vidpath)
%% body_tracker: tracks the body angle of an insect in a magnetic tether
%
%   Fits an ellipse to the 'on' reigon in each frame (insect body). Blurs the
%   image to get rid of head movements and wraps the angle. Can debug by
%   displaying a tracking animation.
%
% Sign convention for angle outputs: [CW = + , CCW = -]
% 0 deg is when the fly hed is pointing upward in the frame
%
%   INPUT:
%       vid             : raw grayscale video data in matrix form
%                         or a path to a video file (.mp4 or .avi)
%       save_path       : path to save calculatd body angles.
%                         if empty, default is 'tracked_body' folder in 'vid' directory
%       mat_var_name    : name of the matlab variable containing the video data
%                         only required if inputing vid_path to .mat file
%       debug_heading   : show debug figure for setting initial heading
%                         0=never, 1=always, 2=if close call (default = 0)
%       flip_vid        : if true, then flip the registered video from left to right
%                         (usually only set to true if a bottom view video)
%       playback   	    : playback rate (show a frame in increments of "playback").
%                         if false, then only show the 1st frame. (default = 1)
%       vidpath         : path of output video. No ouput if false
%                           
%   OUTPUT:
%       angle 	        : normalized & unwrapped body angle [deg]
%       imgstats 	    : structure containing some basic image statistics (orientation, centroid, etc.)
%       initframe       : initial (1st) frame image
%

% Parse inputs
if nargin < 7
    vidpath = [];
    if nargin < 6
        playback = 1;
        if nargin < 5
            flip_vid = false;
            if nargin < 4
                debug_heading = 0;
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
[viddata] = get_vid_data(vid, save_path, mat_var_name, 'tracked_body');

% Flip 1st frame left to right if specified
if flip_vid
    viddata.first_frame = fliplr(viddata.first_frame);
end

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

% Parse playback input
if isempty(playback)
    playback = 1; % default
end

if ~rem(playback,1)==0
    warning('Warning: "playback" rounds to nearest integer')
    playback = round(playback);
end

% Use the inital frame to find the heading
[~, flip] = find_heading(viddata.first_frame, debug_heading);
disp('Heading found.')
pause(0.5)
close all

% Video preprocessing morphological structuring elements
SE_erode = strel('disk', 8, 4); % erosion mask
SE_dilate = strel('disk', 12, 4); % dilation mask

% Set some parameters
offset = 90;  % shift the reference frame so 0 deg is the top vertical axis in video [deg]
dthresh = 160; % threshold for detecting >180 deg flips in ellipse orientation, or angles > 360 deg
shift  = 0; % shift to keep angle wrapped [deg] (dynamic)

% Create display
fig(1) = figure (100); clf
% set(fig, 'Visible', 'off')
fColor = 'k'; % figure and main axis color
aColor = 'w'; % axis line color
set(fig, 'Color', fColor, 'Units', 'inches', 'Position', [2 2 9 7])
fig(1).Position(3:4) = [9 7];
% movegui(fig,'center')
figure (100)
    % Raw image window
    ax(1) = subplot(3, 2, [1,3]); hold on ; cla ; axis image
    H(1) = imshow(viddata.first_frame);

    % Processed video window
    ax(2) = subplot(3, 2, [2,4]); hold on ; cla ; axis image
    H(2) = imshow(viddata.first_frame);

    % Body angle window
    ax(3) = subplot(3,2,5:6); hold on ; cla ; xlim([0 viddata.n_frame])
        xlabel('Frame')
        ylabel('Angle (\circ)', 'Interpreter', 'tex')
        h.raw_angle  = animatedline(ax(3),'Color', 'b', 'LineWidth',1); % debugging
        h.norm_angle = animatedline(ax(3),'Color', 'r', 'LineWidth',1);

set(ax, 'Color', fColor, 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
    'YColor', aColor, 'XColor',aColor)

% Preallocate vectors to store tracked angles
raw_ang = nan(viddata.n_frame,1); % stores raw angles calulated by ellipse fit [deg]
norm_ang = nan(viddata.n_frame,1); % stores normalized/unwrapped angles [deg]

if export
    VID = VideoWriter(vidpath,'MPEG-4');
    VID.Quality = 100;
    VID.FrameRate = 60;
    open(VID)
end

tic
disp('Tracking...')
for n = 1:viddata.n_frame
    % Display frame count every every 100 frames
    if (n==1) || ~mod(n,100) || (n==viddata.n_frame)
        fprintf('%i\n',n)
    end
    
    % Get frame
    if viddata.matflag % from matrix
        frame = viddata.vid(:, :, n);
    else % from video reader 
        frame = im2gray(read(viddata.Vread, n));
    end

    % Invert if backinglighting was detected
    if viddata.invertflag
        frame = imcomplement(frame);
    end

    % Flip video left to right if specified
    if flip_vid
        frame = fliplr(frame);
    end

    % Preprocess
    bnframe = imbinarize(frame); % binarize
    bnframe = imerode(bnframe, SE_erode); % erode
    bnframe = imdilate(bnframe, SE_dilate); % dilate
    bnframe = logical(bnframe); % store bianary frame
    
  	% Fit an ellipse to the frame
    tempstats = regionprops(bnframe,'Centroid','Area','BoundingBox','Orientation', ...
        'MajorAxisLength','MinorAxisLength'); % image reigon properties
    
    % Make sure we only use the largest object in the frame
    [~,sort_area] = sort([tempstats.Area],'descend');
    bodyI = sort_area(1);
    imgstats(n) = tempstats(bodyI); % get stats
    
    % Get body centroid & angle
    centroid = imgstats(n).Centroid; % centroid of image
    L = imgstats(n).MajorAxisLength / 2; % long axis of image
    raw_ang(n) = (imgstats(n).Orientation - offset); % raw angle [deg]
    norm_ang(n) = -(imgstats(n).Orientation - offset + shift); % normalized, unwrapped angle [deg]
 	
    % Check for changes in angle > 180 deg. Correct for the ellipse fit and unwrap angles.
    if n > 1
        dbody = norm_ang(n) - norm_ang(n-1); % change in body angle between frames [deg]
        magd = abs(dbody); % magnitude of change [deg]
        signd = sign(dbody); % direction of change
        if magd > dthresh % 180 or more flip
            flipn = round(magd/180); % how many 180 deg shifts we need
         	shift = -signd*flipn*180; % shift amount [deg]
            norm_ang(n) = norm_ang(n) + shift; % normalized, unwrapped angle [deg]
        end
    elseif n == 1
        % Flip heading by 180 deg if the heading is in the wrong direction
        if flip
            norm_ang(n) = norm_ang(n) + 180;
        end
        
       	% Make start angle positive
        if norm_ang(n) > 180
            norm_ang(n) = norm_ang(n) - 360;
        end
    end

    if playback || (n == 1)
        % Display images
        if (n==1) || (~mod(n,playback)) || (n==viddata.n_frame) % display at playback rate
            % Get approximate location of head
            head = centroid + L*[sind(norm_ang(n)), -cosd(norm_ang(n))];
            heading = [centroid ; head];
            abdomen = centroid - L*[sind(norm_ang(n)), -cosd(norm_ang(n))];
            reverse = [centroid ; abdomen];
            
            % Show images with tracking annotation
            set(fig, 'CurrentAxes', ax(1)) ; cla
                imshow(frame) % frame
                %set(H, 'CData', frame);
    
                h.heading = semi_ellipse(centroid, L, 0.5, 0.90, 180 - norm_ang(n), 'r');
                delete([h.heading{2:4}])
                alpha(h.heading{1},0.3)
                h.heading{1}.LineStyle = 'none';
                
            set(fig, 'CurrentAxes', ax(2)); cla % processed
                imshow(bnframe) % frame

                % Show bounding box
                h.rect = rectangle('Position', imgstats(n).BoundingBox, 'EdgeColor', 'g', ...
                                                                'LineWidth', 1);

                % Show ellipse
                h.ellps = ellipse(centroid, 2*L, 0.5, 0.90, 180 - norm_ang(n), 'r');
                delete([h.ellps{[1,3:6]}])
                h.ellps{2}.Color = 'r';
                h.ellps{2}.LineWidth = 1;
                
                h.heading = semi_ellipse(centroid, L, 0.5, 0.90, 180 - norm_ang(n), 'r');
                delete([h.heading{2:4}])
                alpha(h.heading{1},0.3)
                h.heading{1}.LineStyle = 'none';

                % Show centroid & heading
                plot(heading(:,1), heading(:,2), '-r', 'LineWidth', 1) % centroid
                plot(reverse(:,1), reverse(:,2), '-b', 'LineWidth', 1) % heading
                plot(centroid(1),  centroid(2),  '.g', 'MarkerSize', 20) % reverse heading
        end
        
        % Display angle
        set(fig, 'CurrentAxes', ax(3))
            addpoints(h.norm_angle, n, norm_ang(n))
            
        if n == 1 % get 1st frame
            initframe = getframe(fig);
            initframe = initframe.cdata;
        end
        
        drawnow
        
        if export
          	fig_frame = getframe(fig);
         	writeVideo(VID, fig_frame); 
        end
    end
end
if export
   close(VID) 
end

% Save image data & body angles
angle = norm_ang;
save(viddata.save_data_path, 'angle', 'imgstats', 'initframe', '-v7.3')

toc
end