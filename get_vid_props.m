function [yaw, vid_props, bound_box] = get_vid_props(vid, debug)
%% get_vid_props: get properties of video
%
%   INPUT:
%       vid         :   video matrix
%       debug       :   show comparison for each frame
%
%   OUTPUT:
%       yaw         :   angle of largets object
%       vid_props  	:   video properties
%       bound_box   :  	median bounding box
%

vid = squeeze(vid);
se = strel('disk',7);
for n = 1:size(vid,3)
    A = cell(1);
    A{1} = vid(:,:,n);
    A{end+1} = medfilt2(A{end} ,[3 3]);
    A{end+1} = imbinarize(A{end});
    A{end+1} = bwareaopen(A{end},50);
    A{end+1} = imdilate(A{end}, se);

    image_props = regionprops(A{end},'Centroid','MajorAxisLength','MinorAxisLength',...
                                'BoundingBox','Orientation','Image');

    if length(image_props) > 1
        mjaxes = [image_props.MajorAxisLength];
        [~,idx] = sort(mjaxes, 'descend');
        idx = idx(1);
    else
        idx = 1;
    end                  
    vid_props(n) = image_props(idx);

    if debug
       montage(A)
       pause
    end
end
bound_box = median(cat(1,vid_props.BoundingBox),1);
yaw = [vid_props.Orientation];
yaw(yaw < -40) = yaw(yaw < -40) + 180;
yaw = yaw - 90;

end