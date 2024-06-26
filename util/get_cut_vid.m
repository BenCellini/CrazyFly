function [VID, cent] = get_cut_vid(vid, percent_frame, bb, flex, cut)
% get_head_vid: 
%
%   INPUT:
%       vid     : video matrix
%       P       : what percent of frames to grab (1 < size(P,3))    
%       bb      : bounding box
%       flex    : 1x2 proportion to extend frame by in each direction
%       cut     : proportion to cut frame from bottom
%
%   OUTPUT:
%       VID     : video output structure
%       cent    : median centroid
%       

vid = squeeze(vid);
dim = size(vid);

% Get ROI
if isempty(bb)
    VID.main = vid(:,:,:);
else
    x1 = floor((1 - flex(1)) * bb(1));
    y1 = floor((1 - flex(2)) * bb(2));
    x2 = ceil((1 + flex(1)) * (bb(1) + bb(3)));
    y2 = ceil((1 - cut) * (bb(2) + bb(4)));
    xr = x1:x2;
    yr = y1:y2;
    out_range_x = (xr > dim(2)) | (xr < 1);
    out_range_y = (yr > dim(1)) | (yr < 1);
    xr = xr(~out_range_x);
    yr = yr(~out_range_y);
    VID.main = vid(yr,xr,:);
end

if ~isempty(percent_frame)
	n_use_frame = round(percent_frame*dim(3));
    rand_frame = randperm(dim(3), n_use_frame);
    VID.main = VID.main(:,:,rand_frame);
end
n_frame = size(VID.main, 3);

VID.bw = false(size(VID.main)); % bianarized video
VID.out = false(size(VID.main)); % outline video
se_dilate = strel('disk',8);
se_erode = strel('disk',4);
se_close = strel('disk',1);
for n = 1:n_frame
    A = cell(1);
    A{1} = VID.main(:,:,n);
    A{end+1} = medfilt2(A{end}, [5 5]);
    A{end+1} = imbinarize(A{end});
    A{end}(end,:) = true;
    A{end+1} = imclose(A{end}, se_close);
    A{end+1} = imfill(A{end}, 'holes');
    A{end+1} = bwareaopen(A{end},300);
    
    cent_frame = imdilate( imdilate(A{end}, se_dilate), se_dilate);
    image_props = regionprops(cent_frame,'Centroid','MajorAxisLength','PixelList', 'Area');               
    if length(image_props) > 1
        mjaxes = [image_props.Area];
        [~,idx] = sort(mjaxes, 'descend');
        keepI = idx(1);
        xy_mask = image_props(keepI).PixelList;
        mask = false(dim(1), dim(2));
        for p = 1:size(xy_mask,1)
            mask(xy_mask(p,2), xy_mask(p,1)) = true;
        end
        A{end}(~mask) = false;
    else
        keepI = 1;
    end                  
    cut_props(n) = image_props(keepI);

    % Get outline
    A{end+1} = imerode(A{end}, se_erode);
    A{end+1} = imfill(A{end}, 'holes');
    A{end+1} = bwareaopen(A{end},300);
    A{end+1} = bwmorph(A{end}, 'bridge',3);
    
    VID.bw(:,:,n) = A{end};

    A{end+1} = bwmorph(A{end},'remove');

    VID.out(:,:,n) = A{end};
end
cent = median(cat(1,cut_props.Centroid),1);
end