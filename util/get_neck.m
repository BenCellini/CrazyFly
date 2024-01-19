function [pivot, R_head, R_body, body_angle] = get_neck(vid, vid_bw, cent)
%% get_neck: estimate neck point
%
%   INPUT:
%       vid         :   video matrix
%       vid_bw    	:   bianarized video
%       cent        :   centroid
%
%   OUTPUT:
%    	pivot       :   neck pivot point
%       R_head     	:   head radius
%       R_body     	:   body radius
%     	body_angle 	:   body angle [°]
%

dim = size(vid);
Fc_n = 0.2; % normalized filtering cut-off for smoothing fly outline
for n = 1:dim(3)
    frame = vid(:,:,n); % frame for analysis

    % Get left and right segments of frame
    left = frame;
    left(:,round(cent(1)):dim(2)) = false;
    right = frame;
    right(:,1:round(cent(1))) = false;

    % Get the farthest 'on' pixel from center line for left and right images
    xx = zeros(dim(1),2);
    xR = zeros(dim(1),1);
    xL = zeros(dim(1),1);
    for r = 1:dim(1)
        all_x = find(frame(r,:));
        if length(all_x) > 1
            xx(r,:) = [all_x(1) all_x(end)];
        end
        x_right = find(right(r,:), 1, 'last');
        x_left = find(left(r,:), 1, 'first');
        if ~isempty(x_right)
            xR(r) = find(right(r,:), 1, 'last');
        end
        if ~isempty(x_left)
            xL(r) = find(fliplr(left(r,:)), 1, 'last');
        end
    end

    % Filter curves
    [b, a] = butter(3, Fc_n, 'low');
    xR_filt = filtfilt(b, a, xR);
    xL_filt = filtfilt(b, a, xL);

    % Get neck properties in frame
    if n > 1
        [nR(n)] = get_neck_point(xR_filt, right, nR, false);
        [nL(n)] = get_neck_point(xL_filt, left, nL, false);
    else
        [nR(n)] = get_neck_point(xR_filt, right, [], false);
        [nL(n)] = get_neck_point(xL_filt, left, [], false);
    end
end

% Estimate the y-axis pivot & radius of the neck about yaw
pivot(2) = 0.96*nanmedian(cat(1, [nL.peak]', [nR.peak]'));
startI_med = median(cat(1, [nL.startI]', [nR.startI]'));
endI_med = median(cat(1, [nL.endI]', [nR.endI]'));
neck_y_range = round( (pivot(2) - startI_med) : pivot(2));
body_y_range = round(pivot(2)+10:endI_med);
R_head = 1.1*(pivot(2) - startI_med);
R_body = (endI_med - pivot(2));

if any(neck_y_range <= 0)
    warning('neck detetcion unsuccessful')
    neck_y_range(neck_y_range <= 0) = [];
end

if any(body_y_range <= 0)
    warning('neck detetcion unsuccessful')
    body_y_range(body_y_range <= 0) = [];
end

% Get the x-axis pivot, body centroid, & body angle
neck_x_all = nan(dim(3),1);
body_cent = nan(dim(3),2);
for n = 1:dim(3)
   frame = vid_bw(neck_y_range,:,n);
   [~,neck_x] = find(frame);
   neck_x_all(n) = mean(neck_x);
   
   frame_body = vid_bw(:,:,n);
   frame_body(1:body_y_range(2),:) = false;
   [body_y,body_x] = find(frame_body);
   body_cent(n,:) = median([body_x body_y]);
end
pivot(1) = median(neck_x_all);
body_cent = median(body_cent,1);
body_axis = body_cent - pivot;
body_angle = -rad2deg(atan2(body_axis(1), body_axis(2)));

end

function [neck] = get_neck_point(edge, I, side, debug)
%% get_neck_point:
%   INPUTS:
%       edge    : far edge for each row
%       I       : left or right image
%       side    : structure with previous neck points
%       debug   : show plots
%
%   OUTPUTS:
%       neck    : neck structure
%

neck.I = I;
neck.full_edge = edge;
startI = round(1.1*find(any(I,2),1,'first'));
endI = round(0.7*find(any(I,2),1,'last'));
cut_edge = nan(size(edge));
cut_edge(startI:endI) = edge(startI:endI);

neck.startI = startI;
neck.endI = endI;

[pks,locs,w,p] = findpeaks(-cut_edge, 'MinPeakDistance', 10, 'MinPeakWidth', 3, ...
                                'MinPeakProminence', 10, 'SortStr', 'descend');
idx = (1:length(edge))';
if isempty(locs) % didn't find the neck on first try
    %warning('Small peak')
    neck.flag = true;
    if ~isempty(side) % use the last neck point
        good_peak = ~[side.flag]';
        pks = -[side.dx];
        pks = median(pks(good_peak));
        locs = [side.peak];
        locs = median(locs(good_peak));
        w = [side.width];
        w = median(w(good_peak));
    else % lower peak detection parameter thresholds
        [pks,locs,w,p] = findpeaks(-edge, 'MinPeakDistance', 10, ...
                'MinPeakWidth', 3, 'MinPeakProminence', 5);
    end

    if isempty(locs) % if we still can't find the peak then set to NaN
        warning('Can''t find angle')
        neck.flag = nan;
        neck.peak = nan;
        neck.dx = nan;
        neck.width = nan;
        return
    else
        locs = locs(end);
        pks = pks(end);
    end
else
    neck.flag = false;
end

neck.peak = round(locs(1));
neck.dx = -pks(1);
neck.width = w(1);

if debug
    figure (109); cla
    ax(1) = subplot(1,1,1); hold on; axis image
        imshow(255*I, 'InitialMagnification', 400) ; hold on
        plot(edge,idx,'m', 'LineWidth', 1)
        plot(cut_edge,idx,'r', 'LineWidth', 1)
        plot(neck.dx, neck.peak, '.g', 'MarkerSize', 15)
    set(ax, 'Box', 'on')
end
end
        