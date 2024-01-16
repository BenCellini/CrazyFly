function [stable_iso_vid, cent, vid_props] = stable_head(vid, yaw, pivot)
%% stable_head: isolate head in yaw stabilized video
%
%   INPUT:
%       vid         :   video matrix
%       yaw         :   yaw angle in each frame
%       pivot       :   head pivot point
%
%   OUTPUT:
%       stable_vid 	:   stabilized video
%       cent        :   median centroid
%       vid_props   :   stable iso video properties
%

vid = squeeze(vid);
dim = size(vid);

[b,a] = butter(3, 0.5, 'low');
yaw_filt = yaw;
yaw_filt = hampel(1:dim(3), yaw_filt, 10, 3, 'Adaptive', 0.1);
yaw_filt = filtfilt(b, a, yaw_filt);

% Stabilize the frame in head reference
stable_vid = uint8(false(size(vid)));
for n = 1:dim(3)
    stable_vid(:,:,n) = imrotate(vid(:,:,n), ...
        yaw_filt(n), 'bilinear', 'crop');
end

% Find where the head starts in every frame (from top and sides)
head_top = nan(dim(3),1);
head_lr = nan(dim(3),2);
for n = 1:dim(3)
    frame = vid(:,:,n);
    frame_bw = imbinarize(frame);
    frame_bw = bwareaopen(frame_bw, 5);
    head_top(n) = find( max(frame_bw,[],2) , 1, 'first');
    
    if n < pivot(2)
        head_lr(n,1) = find( any(frame_bw, 1), 1, 'first');
        head_lr(n,2) = find( any(frame_bw, 1), 1, 'last');
    end
end
head_top = prctile(head_top,10);
head_tb = round(head_top:pivot(2));
head_lr = round(prctile(head_lr,50,1));
head_lr = head_lr(1):head_lr(2);

h = length(head_tb);
w = length(head_lr);

% Isolate head in each frame
stable_iso_vid = uint8( zeros( h, w) );
for n = 1:dim(3)
    stable_iso_vid(:,:,n) = stable_vid(head_tb,head_lr,n);

    bw_iso = imbinarize(stable_iso_vid(:,:,n));
    image_props = regionprops(bw_iso,'Centroid','BoundingBox','MajorAxisLength');

    if length(image_props) > 1
        mjaxes = [image_props.MajorAxisLength];
        [~,idx] = sort(mjaxes, 'descend');
        idx = idx(1);
    else
        idx = 1;
    end  

    vid_props(n) = image_props(idx);
end
cent = median(cat(1, vid_props.Centroid),1);

end