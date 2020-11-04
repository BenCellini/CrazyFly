function [roll, roll_idx, eye] = track_head_roll(vid, yaw, pivot, eye_ratio, roll_cal, debug)
% get_roll: find roll angle of head in each frame

%  roll_cal: roll calibration factor
%
%   INPUT:
%       vid         :   video matrix
%       yaw         :   yaw angle in each frame
%       pivot       :	head pivot point
%       eye_ratio   :   ratio of frame height around centroid to get eye intensity
%       roll_cal    :   roll calibration factor
%
%   OUTPUT:
%       stable_vid 	:   stabilized video
%       cent        :   median centroid
%       vid_props   :   stable iso video properties
%

vid = squeeze(vid);
dim = size(vid);

% Get stabilized head window
[stable.vid, stable.cent] = stable_head(vid, yaw, pivot);
res = 5;
dim_stable = size(stable.vid);

% Find window to measure eye intensity
eye.xspan = 1:dim_stable(2);
eye.xinterp = 1:1/res:dim_stable(2);
eye.cent = stable.cent;
eye.span = round(eye_ratio * dim_stable(1) / 2);
eye.window = round((eye.cent(2) - eye.span : eye.cent(2) + eye.span)');
eye.intensity = cell(dim_stable(3),1);

% Find roll angles
[b,a] = butter(3, 0.3, 'low');
for n = 1:dim(3)
    frame = stable.vid(:,:,n); % frame for roll analysis

    % Get eye intensities in window
    eye.intensity{n} = double(frame(eye.window,:)); % intensity in eye range
    int_mean = filtfilt(b, a, mean(eye.intensity{n},1));
    int_mean_interp = interp1(eye.xspan, int_mean, eye.xinterp, 'pchip');
    eye.intensity_mean(n,:) = int_mean_interp; % mean intensity
    dx = diff(eye.intensity_mean(n,:)) / mean(diff(eye.xinterp));
    dx = [dx(1) , dx];
    dx = filtfilt(b, a, dx);
    eye.intensity_mean_diff(n,:) = dx; % derivative of mean intensity

    % Find outer spikes to get eye starts (outer)
    [~,locs,w] = findpeaks(eye.intensity_mean_diff(n,:), 'MinPeakHeight', 2, ...
                    'MinPeakWidth', res*2, 'MinPeakProminence', 0.5);
    left_outer_I = (locs(1) - 0*w(1)/1.5);
    left_outer_x = eye.xinterp(left_outer_I);

    [~,locs,w] = findpeaks(-eye.intensity_mean_diff(n,:), 'MinPeakHeight', 2, ...
                    'MinPeakWidth', res*2, 'MinPeakProminence', 0.5);
    right_outer_I = (locs(end) + 0*w(end)/1.5);
    right_outer_x = eye.xinterp(right_outer_I);

    % Find peak intensities
    [pks,locs,w,p] = findpeaks(eye.intensity_mean(n,:), 'MinPeakHeight', 100, ...
                    'MinPeakDistance', res*5, 'MinPeakWidth', res*2, 'MinPeakProminence', 5);

    % Get left and right eye peaks (1st & last)
    eye.peak_loc(n,:) = locs([1,end]); % pixel in horizontal plane
    eye.peak_x(n,:) = eye.xinterp(locs([1,end])); % interp pixel in horizontal plane
    eye.peak_int(n,:) = pks([1,end]); % intensity value (close too 255)

    % Get change in intensity between peaks
    flex = round(0.05 * dim_stable(2));
    eye.eye_range{n} = eye.peak_loc(n,1)-flex : eye.peak_loc(n,2)+flex; % between peaks
    mid = nan(size(dx));
    mid(eye.eye_range{n}) = eye.intensity_mean_diff(n,eye.eye_range{n});
    eye.eye_mid{n} = mid; % intensity between peaks

    [~,locs,w] = findpeaks(-eye.eye_mid{n}, 'MinPeakHeight', 1, ...
                    'MinPeakWidth', res*2, 'MinPeakProminence', 3);
    if ~isempty(locs)
        left_inner = (locs(1) + 0*w(1)/1.5);
        eye.left_inner_int(n) = eye.intensity_mean(n,round(left_inner));
        eye.left(n,:) = [left_outer_x, eye.xinterp(left_inner)];
    else
        eye.left_inner_int(n) = 0;
        eye.left(n,:) = [0, 0];
    end

    [~,locs,w] = findpeaks(eye.eye_mid{n},'MinPeakHeight', 1, ...
                    'MinPeakWidth', res*2, 'MinPeakProminence', 3);
    if ~isempty(locs)
        right_inner = (locs(end) - 0*w(end)/1.5);
        eye.right_inner_int(n) = eye.intensity_mean(n,round(right_inner));
        eye.right(n,:) = [eye.xinterp(right_inner) , right_outer_x];
    else
        eye.right_inner_int(n) = 0;
        eye.right(n,:) = [0 , 0];
    end

    if debug
        figure (111) ; clf
        subplot(2,1,1)
            imshow(frame) ; hold on
            x = [1 1 dim_stable(2) dim_stable(2)];
            y = [eye.window(1) eye.window(end) eye.window(end) eye.window(1)];
            patch(x, y, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
            plot(eye.left(n,:), [eye.cent(2) eye.cent(2)],  '-r', 'LineWidth', 3)
            plot(eye.right(n,:), [eye.cent(2) eye.cent(2)],  '-r', 'LineWidth', 3)
            
        subplot(2,1,2) ; hold on ; axis tight ; ylim([-50 260])
            plot(eye.xinterp, eye.intensity_mean(n,:), 'k', 'LineWidth', 1)
            plot(eye.xinterp([1 end]), [0 0], '--k', 'LineWidth', 0.5)
            plot(eye.peak_x(n,:), eye.peak_int(n,:), '.r', 'MarkerSize', 15)
            plot(eye.xinterp, eye.intensity_mean_diff(n,:), 'g', 'LineWidth', 1)
            plot(eye.xinterp, eye.eye_mid{n}, 'r', 'LineWidth', 1.5)
            plot(eye.left(n,:), [eye.left_inner_int(n) eye.left_inner_int(n)], ...
                '.-m', 'MarkerSize', 15)
            plot(eye.right(n,:), [eye.right_inner_int(n) eye.right_inner_int(n)], ...
                '.-c', 'MarkerSize', 15)
            plot([eye.left(n,1) eye.left(n,1)], [0 eye.left_inner_int(n)], 'm')
            plot([eye.left(n,2) eye.left(n,2)], [0 eye.left_inner_int(n)], 'm')
            plot([eye.right(n,2) eye.right(n,2)], [0 eye.right_inner_int(n)], 'c')
            plot([eye.right(n,1) eye.right(n,1)], [0 eye.right_inner_int(n)], 'c')
    end
end
% Find eye widths and roll index
eye.left_width  = diff(eye.left, 1, 2);
eye.right_width = diff(eye.right, 1, 2);
roll_idx = (eye.left_width - eye.right_width) ./ ...
    (eye.left_width + eye.right_width);
roll = roll_idx .* roll_cal;

end   