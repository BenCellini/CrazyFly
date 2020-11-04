function [mask] = ui_arc_mask(vidData)
%%
clearvars -except vidData regvid t_v hAngles
clc
close all
vid = squeeze(vidData);
dim = size(vid);
frame = vid(:,:,100);

%% Bouding box
[body_yaw, vid_props, bb] = get_vid_props(vid, false);

%% Head vid
[VID,cent] = get_cut_vid(vid, [], [0.1 0.1], 0.3);

%% Neck joint
[pivot,R_head] = get_neck(VID.out, VID.bw, cent);

%% Mask
% "H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\30\Vid\
% Fly_10_Trial_21_Vel_120_SpatFreq_30.mat"
frame = vid(:,:,100);

global mask
% mask = make_mask(pivot, median(body_yaw), [R_head 1.75*R_head], [40 40], VID.main(:,:,1));
mask = make_mask(pivot, 0, [R_head 1.75*R_head], [40 40], VID.main(:,:,1));
% uiwait(mask.fig.main)

disp('Done')

%%
clc
figure (1)
clf
npts = 120;
norm = 2;
yaw = nan(size(vid,3),1);
M = nan(size(vid,3),2);
for n = 1:size(vid,3)
    %frame = vid(:,:,n);
    frame = VID.main(:,:,n);
    [angle,m,pts,k] = tracktip(frame, mask.area_points, ...
        mask.move_points.rot, norm, npts, 'clust');
    M(n,:) = m';
    angle = angle - 270;
    rot = mask.move_points.rot;
    tip =  rot + 45*[sind(angle) -cosd(angle)];
    
    yaw(n) = angle - mask.global;
    
	imshow(frame) ; hold on ; title(angle)
    plot(pts(k==1,1),pts(k==1,2), 'r.', 'MarkerSize', 5)
    plot(pts(k==2,1),pts(k==2,2), 'b.', 'MarkerSize', 5)
    plot([rot(1) tip(1)], [rot(2) tip(2)], 'Color', mask.color, 'LineWidth', 2)
    
	pause(0.001)
    %pause
    cla
end
disp('Done')

%%
[stable_iso_vid, cent] = stable_head(vid, yaw, pivot);

end

function mask = make_mask_old(rot, global_ang, R, span, frame)
    global H mask
    
	% Defaults
    mask.lw = 0.1;
    mask.color = [0 1 1]; % mask color
    mask.sym = true;
    
    % Main window
    mask.fig.main = figure('Name', 'Flyalyzer', 'Color', [0.3 0.3 0.3]);
    imshow(frame)
    
    fig_w = mask.fig.main.Position(3);
	fig_h = mask.fig.main.Position(4);
    
    % UI Controls
    H.sym = uicontrol('style','checkbox','units','pixels', 'Value', mask.sym, ...
                    'position',[10,0.05*fig_h,70,15],'string','symmetric');
    set(H.sym, 'callback', @(src, event) sym_call(src, event, H));

    H.zero = uicontrol('style','edit','units','pixels', ...
                    'position',[10,0.1*fig_h,50,15],'string','0');
    set(H.zero, 'callback', @(src, event) zero_call(src, event, H));
    
    H.done = uicontrol('style','pushbutton','units','pixels', 'Value', false, ...
                    'position',[10,0.01*fig_h,50,15],'string','Done');
    set(H.done, 'callback', @(src, event) done_call(src, event, H));
    
    % Get frame to set mask
    mask.image = frame;
    mask.image_sz = size(mask.image);
    
    % Intiial mask settings
    mask.global = global_ang; % body angle (mask reference frame) [°]
    mask.radius.inner = R(1); % inner mask radius
    mask.radius.outer = R(2); % outer mask radius
    mask.radius.axis = mask.radius.outer; % axis radius
    mask.span = span; % angular span of mask [+,-] [°]
    
    % Define rotation point
    hold on
    mask.move.rot = drawpoint('Position', rot, 'Color', mask.color);
    mask.mrk_sz = 0.7 * mask.move.rot.MarkerSize;
    mask.move.rot.MarkerSize = mask.mrk_sz;
    
    % Define mask angles in global frame
    mask.left_angle = mask.global - mask.span(1); % left angle in global [°]
    mask.right_angle = mask.global + mask.span(2); % left angle in global [°]
    mask.angles = (mask.left_angle:1:mask.right_angle)'; % angles in between edges of mask
    
    % Define arc coordinates
    inner_cent_pos = mask.move.rot.Position + ...
        mask.radius.inner*[sind(mask.global) , -cosd(mask.global)];
    outer_cent_pos = mask.move.rot.Position + ...
        mask.radius.outer*[sind(mask.global) , -cosd(mask.global)];
    outer_left_pos = mask.move.rot.Position + ...
        mask.radius.outer*[sind(mask.left_angle) , -cosd(mask.left_angle)];
    outer_right_pos = mask.move.rot.Position + ...
        mask.radius.outer*[sind(mask.right_angle) , -cosd(mask.right_angle)];
    axis_pos = mask.move.rot.Position - ...
        4*mask.radius.axis*[sind(mask.global) , -cosd(mask.global)];
    
    % Make the movable mask ROI points
    mask.move.inner_C = drawpoint('Position', inner_cent_pos, ...
        'Color', mask.color, 'MarkerSize', mask.mrk_sz);    
    mask.move.outer_C = drawpoint('Position', outer_cent_pos, ...
        'Color', mask.color, 'MarkerSize', mask.mrk_sz);
    mask.move.outer_L = drawpoint('Position', outer_left_pos, ...
        'Color', mask.color, 'MarkerSize', mask.mrk_sz);
    mask.move.outer_R = drawpoint('Position', outer_right_pos, ...
        'Color', mask.color, 'MarkerSize', mask.mrk_sz);
    mask.move.axis = drawpoint('Position', axis_pos, ...
        'Color', mask.color, 'MarkerSize', mask.mrk_sz);
    
    % Listeners to move mask
    addlistener(mask.move.rot, 'MovingROI', @rot_point);
    addlistener(mask.move.rot, 'ROIMoved', @rot_point);
    addlistener(mask.move.inner_C, 'MovingROI', @inner_C);
    addlistener(mask.move.inner_C, 'ROIMoved', @inner_C);
    addlistener(mask.move.outer_C, 'MovingROI', @outer_C);
    addlistener(mask.move.outer_C, 'ROIMoved', @outer_C);
    addlistener(mask.move.outer_L, 'MovingROI', @outer_left);
    addlistener(mask.move.outer_L, 'ROIMoved', @outer_left);
    addlistener(mask.move.outer_R, 'MovingROI', @outer_right);
    addlistener(mask.move.outer_R, 'ROIMoved', @outer_right);
    addlistener(mask.move.axis, 'MovingROI', @axis);
    addlistener(mask.move.axis, 'ROIMoved', @axis);
    
    % Draw mask
    mask.patch = [];
    mask.axis = [];
    mask.axis_reverse = [];
    mask = draw_mask(mask);
    mask = get_ROI(mask);
end

function mask = draw_mask(mask)
    % Define new mask angles
    mask.span(mask.span < 0) = mask.span(mask.span < 0) + 360;
    mask.angles = (mask.global + ( -mask.span(1):mask.span(2) ))';
    mask.left_angle = mask.angles(1);
    mask.right_angle = mask.angles(end);

    % Create the new arc mask point map
    mask.inner.pts = mask.radius.inner*[sind(mask.angles) , -cosd(mask.angles)];
    mask.outer.pts = mask.radius.outer*[sind(mask.angles) , -cosd(mask.angles)];
    mask.points = mask.move.rot.Position + ...
        [mask.inner.pts ; flipud(mask.outer.pts) ; mask.inner.pts(1,:)];
     
	% Constrain points
    mask.move.inner_C.Position = mask.move.rot.Position + ...
        mask.radius.inner*[sind(mask.global) , -cosd(mask.global)];
    mask.move.outer_C.Position = mask.move.rot.Position + ...
        mask.radius.outer*[sind(mask.global) , -cosd(mask.global)];
    mask.move.outer_L.Position = mask.move.rot.Position + ...
        mask.radius.outer*[sind(mask.left_angle) , -cosd(mask.left_angle)];
    mask.move.outer_R.Position = mask.move.rot.Position + ...
        mask.radius.outer*[sind(mask.right_angle) , -cosd(mask.right_angle)];
    mask.move.axis.Position = mask.move.rot.Position - ...
        mask.radius.axis*[sind(mask.global) , -cosd(mask.global)];
    
    % Annotate the mask ROI reigon
    delete(mask.patch)
    delete([mask.axis mask.axis_reverse])
    mask.patch = patch(mask.points(:,1), mask.points(:,2), ...
        mask.color, 'FaceAlpha', 0.2, 'EdgeColor', mask.color, 'LineWidth', 1);
    mask.axis = plot([mask.move.rot.Position(1) mask.move.outer_C.Position(1)], ...
        [mask.move.rot.Position(2) mask.move.outer_C.Position(2)], ...
         '--', 'Color', mask.color', 'LineWidth', 0.25);
    mask.axis_reverse = plot([mask.move.axis.Position(1) mask.move.rot.Position(1)], ...
        [mask.move.axis.Position(2) mask.move.rot.Position(2)], ...
         '-', 'Color', mask.color', 'LineWidth', 0.25);
     
    % Store the points
    mask.move_points = structfun(@(x) x.Position, mask.move, 'UniformOutput', false);
     
	% Update global heading display
	global H
	H.zero.String = string(mask.global);
end

function mask = get_ROI(mask, frame, showplot)
    if nargin < 3
       showplot= false; 
       if nargin < 2
          frame = mask.image; 
       end
    end

  	% Get mask reigon
	mask.area_points = poly2mask(mask.points(:,1), mask.points(:,2), ...
        mask.image_sz(1), mask.image_sz(2));
    mask.ROI = round([min(mask.points(:,1)) max(mask.points(:,1)), ... 
        min(mask.points(:,2)) max(mask.points(:,2))]);
    mask.ROI(mask.ROI < 0) = 1;
    if mask.ROI(2) > mask.image_sz(2)
        mask.ROI(2) = mask.image_sz(2);
    end
    if mask.ROI(4) > mask.image_sz(1)
        mask.ROI(4) = mask.image_sz(1);
    end
    
    % Get mask image
    if nargin == 1
        mask.frame = mask.image;
    elseif nargin==2
        mask.frame = frame;
    end        
    %mask.frame_bw = imbinarize(mask.frame);
    mask.frame = frame;
    mask.bw_thresh = graythresh(mask.frame(mask.area_points));
    mask.frame(~mask.area_points) = 0;
    %mask.frame_bw(~mask.area_points) = 0;
    mask.y = round(mask.ROI(3):mask.ROI(4));
    mask.x = round(mask.ROI(1):mask.ROI(2));
    %mask.xy = mask.y
    mask.frame = mask.frame(mask.y, mask.x);
    %mask.frame_bw = mask.frame_bw(mask.y, mask.x);
    mask.frame_bw = imbinarize(mask.frame, mask.bw_thresh);
    
    if showplot
        fig = figure(2);
        set(fig, 'Name','ROI')
        imshow(mask.frame_bw)
    end
end

function rot_point(src,evt)
    global mask
    evname = evt.EventName;
    switch(evname)
        case{'MovingROI'}
            shift = evt.CurrentPosition - evt.PreviousPosition;
            mask.move.inner_C.Position = mask.move.inner_C.Position + shift;
            mask.move.outer_C.Position = mask.move.outer_C.Position + shift;
            mask.move.outer_L.Position = mask.move.outer_L.Position + shift;
            mask.move.outer_R.Position = mask.move.outer_R.Position + shift;
            mask = draw_mask(mask);
        case{'ROIMoved'}
            mask = get_ROI(mask);
    end
end

function axis(src,evt)
    global mask
    evname = evt.EventName;
    switch(evname)
        case{'MovingROI'}
            shift = evt.CurrentPosition - mask.move.rot.Position;
            mask.radius.axis = sqrt(shift(1)^(2) + shift(2)^(2));
            mask.global = rad2deg(atan2(shift(2), shift(1))) - 90;
            mask = draw_mask(mask);
        case{'ROIMoved'}
            mask = get_ROI(mask);
    end
end

function outer_C(src,evt)
    global mask
    evname = evt.EventName;
    switch(evname)
        case{'MovingROI'}
            shift = evt.CurrentPosition - mask.move.rot.Position;
            mask.radius.outer = sqrt(shift(1)^(2) + shift(2)^(2));
            %mask.global = rad2deg(atan2(shift(2), shift(1))) + 90;
            mask = draw_mask(mask);
        case{'ROIMoved'}
            mask = get_ROI(mask);
    end
end

function outer_left(src,evt)
    global mask
    evname = evt.EventName;
    switch(evname)
        case{'MovingROI'}
            shift = mask.move.rot.Position - evt.CurrentPosition;
            mask.span(1) = rad2deg(atan2(shift(1), shift(2))) + mask.global;
           	if mask.sym
                mask.span(2) = mask.span(1);
            end
            
            shift = evt.CurrentPosition - mask.move.rot.Position;
            mask.radius.outer = sqrt(shift(1)^(2) + shift(2)^(2));           
            mask = draw_mask(mask);
        case{'ROIMoved'}
            mask = get_ROI(mask);
    end
end

function outer_right(src,evt)
    global mask
    evname = evt.EventName;
    switch(evname)
        case{'MovingROI'}
            shift = mask.move.rot.Position - evt.CurrentPosition;
            mask.span(2) = rad2deg(atan2(-shift(1), shift(2))) - mask.global;
           	if mask.sym
                mask.span(1) = mask.span(2);
            end
            
            shift = evt.CurrentPosition - mask.move.rot.Position;
            mask.radius.outer = sqrt(shift(1)^(2) + shift(2)^(2));
            mask = draw_mask(mask);
        case{'ROIMoved'}
            mask = get_ROI(mask);
    end
end

function inner_C(src,evt)
    global mask
    evname = evt.EventName;
    switch(evname)
        case{'MovingROI'}
            shift = evt.CurrentPosition - mask.move.rot.Position;
            mask.radius.inner = sqrt(shift(1)^(2) + shift(2)^(2));
            mask = draw_mask(mask);
        case{'ROIMoved'}
            mask = get_ROI(mask);
    end
end

% Symmetric checkbox callback
function sym_call(src, event, H)
    global mask
    mask.sym = H.sym.Value;
end

% Zero global callback
function zero_call(src, event, H)
    global mask
    mask.global = str2double(H.zero.String);
    mask = draw_mask(mask);
end

% Done/exit callback
function done_call(src, event, H)
    global mask
    close(mask.fig.main)
end

function clear_mask(mask)
    delete(mask.move.rot)
    delete(mask.move.inner_C)
    delete(mask.move.outer_C)
    delete(mask.move.outer_L)
    delete(mask.move.outer_R)
    delete(mask.patch)
end