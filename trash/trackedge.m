function [edge_ang, tip] = trackedge(img, mask, rot, norm, showplot)
%% get_vid_props: get properties of video
%
%   INPUT:
%       img   	:   image to track tip
%       mask    :   mask logical matrix
%       rot    	:   rotation point
%       norm   	:   normalize contrast based on entire image or just ROI
%       mode  	:   median or k-means
%
%   OUTPUT:
%       angle 	: measured angle [°]
%       tip     : tip location [x,y]
%

if nargin < 5
    showplot = true;
    if nargin < 4
        norm = true;
    end
end

% Convert to greyscale if needed
rot = double(rot);
if size(img,3) > 1
    img = rgb2gray(img);
end
dim = size(img); % size of image

if norm == 2 % normalize tracking ROI only
    img = imadjust(img, stretchlim(img(mask)));
elseif norm == 3 % normalize whole image
   img = imadjust(img);
end

% Get all points in the mask & convert to polar coordinates
[y,x] = find(mask);
I = img;
I(~mask) = 0;
xx = x - rot(1);
yy = -(y - rot(2));
[theta, rho] = cart2pol(xx,yy);
theta = rad2deg(theta);

% Creat radius & angle vectors covering entire mask areea
mask_radii_span = min(rho):0.5:max(rho);
ang_res = 0.5;
mask_ang_span = min(theta):ang_res: max(theta);
if range(mask_ang_span) > 180 % make sure we get the actual mask & not everything outside of it
    sense = sign(min(theta));
    if sense < 0
        mask_ang_span = [min(theta):-ang_res:0 , ang_res: max(theta)];
    end
end

% Get the intensity map in the mask based on the radii & angle of each pixel
n_r = length(mask_radii_span);
n_ang = length(mask_ang_span);
int_map = nan(n_ang,1);
% imshow(I)
% hold on
% scatter(rot(1), rot(2), 50, 'g', 'Marker', '.')
for n = 1:n_ang
    for r = 1:n_r
        pixel_loc = round( rot + (mask_radii_span(r)*.... % location of pixel
            [cosd(mask_ang_span(n)) -sind(mask_ang_span(n))]));
        
        % Make sure pixel is on image
        if pixel_loc(1) > dim(2)
           continue
        elseif pixel_loc(1) < 1
           continue
        end
        if pixel_loc(2) > dim(1)
           continue
        elseif pixel_loc(2) < 1
           continue
        end
        
        pI = img(pixel_loc(2),pixel_loc(1)); % intensity value
        int_map(n,r) = pI; % add value to map
        
%         if mask_radii_span(r) == mask_radii_span(end)
%             plot([rot(1) pixel_loc(1)], [rot(2) pixel_loc(2)], 'g')
%         end
    end
end

% Mean intensity for each angle
int_med = mean(int_map, 2);
int_med_filt = medfilt1(int_med, 5);

% Change in intensity for each angle
[b,a] = butter(5, 0.05);
int_med_diff = central_diff(int_med_filt);
int_med_diff_filt = filtfilt(b, a, int_med_diff);

% Find peaks in intensity gradient
[pks,locs,w,p] = findpeaks(-int_med_diff_filt, ...
    'MinPeakHeight', 0.5, 'MinPeakProminence', 0.2, ...
    'MinPeakWidth', 5, 'MinPeakDistance', 5);

% Pull out the gradient corresponding to the edge
edge_loc = locs(end);
edge_pk = int_med_diff_filt(edge_loc);

% Find the true edge angle
edge_range = int_med_diff_filt;
edge_range(1:edge_loc) = nan;
end_percent = 0.9;
edgeI = find(edge_range > end_percent*edge_pk, 1, 'first');
edge_ang = mask_ang_span(edgeI);

% Tip location
tip = rot + max(rho) * [cosd(edge_ang) -sind(edge_ang)];

if showplot
    fig = figure(203);
    set(fig, 'Name', 'EgdeDetector', 'Color', 'w', 'Units', 'inches', 'Position', [3 1 14 3])
    ax(1) = subplot(1,3,1); hold on ; cla
        imshow(I)
        axis on ; axis image
        plot(rot(1), rot(2), '*m', 'MarkerSize',  8)
        plot([rot(1) tip(1)], [rot(2) tip(2)], 'r', 'LineWidth', 1)
    ax(2) = subplot(1,3,2); hold on ; cla
        [X,Y] = meshgrid(mask_radii_span, mask_ang_span);
        surf(X, Y, int_map, 'EdgeColor', 'none')
        view(2)
        yline(edge_ang, 'r--', 'LineWidth', 1.5)
        xlim([min(mask_radii_span) max(mask_radii_span)])
        ylim([min(mask_ang_span) max(mask_ang_span)])
        xlabel('Radius (pixels)')
        ylabel('Angle (°)')
        axis on
        colormap(ax(2), parula)
    ax(3) = subplot(1,3,3); hold on ; cla
        yline(0, '--', 'Color', [0.5 0.5 0.5])
        plot(mask_ang_span, int_med_diff, 'k', 'LineWidth', 1)
        plot(mask_ang_span, int_med_diff_filt, 'b', 'LineWidth', 1)
        plot(mask_ang_span, edge_range, 'm', 'LineWidth', 1)
        plot(mask_ang_span(locs), -pks, '.r', 'MarkerSize', 15)
        plot(mask_ang_span(edge_loc), edge_pk, 'og', 'MarkerSize', 15)
        xline(edge_ang, 'r--', 'LineWidth', 1.5)
        xlim([min(mask_ang_span) max(mask_ang_span)])
        xlabel('Angle (°)')
        ylabel('Pixel gradient')

    set(ax, 'Color', 'none')
end

end