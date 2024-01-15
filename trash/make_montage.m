function [] = make_montage(vids, labels, outpath)
%% make_montage: makes monatge of two videos of the same size
%
%   INPUT:
%       vids        : cell array of 3D matrix data, one video per cell
%       labels      : title label for each video in cell array
%       outpath     : output video path
%
%   OUTPUT:
%       -
%

n_vid = length(vids);
vids = cellfun(@(x) squeeze(x), vids, 'UniformOutput', false);
norm_size = size(vids{1});

vids{1} = fliplr(vids{1});
for n = 1:n_vid
    for f = 1:norm_size(3)
        vids{n}(:,:,f) = imadjust(vids{n}(:,:,f));
    end
end

for n = 2:n_vid
    sz = size(vids{n});
    padsize = [norm_size(1) - sz(1), norm_size(2) - sz(2)] ./ 2;
    vids{n} = padarray(vids{n}, fix(padsize), 'pre');
    vids{n} = padarray(vids{n}, ceil(padsize), 'post');
end

vid_montage = cat(2, vids{:});

VID = VideoWriter(outpath,'MPEG-4');
VID.FrameRate = 100;
open(VID)

fig = figure;
set(fig, 'Color', 'k', 'Units', 'inches', 'Position', [2 2 6 4])
H = imshow(vid_montage(:,:,1));
for f = 1:norm_size(3)
    if (f == 1) || (f == norm_size(3)) || ~mod(f,10)
        disp(f)
    end
    frame = vid_montage(:,:,f);
    for n = 1:n_vid
        frame = insertText(frame, [(0.5 + (n-1))*norm_size(2) , 25], labels{n}, ...
            'FontSize', 28, 'TextColor', 'white', 'BoxColor', 'red', ...
            'AnchorPoint', 'Center');
    end
    
    set(H, 'CData', frame);
    fig_frame = getframe(fig);
%     fig_frame.cdata = fig_frame.cdata(...
%         round(0.1*norm_size(2)):end-round(0.1*norm_size(2)), ...
%         round(0.01*norm_size(1)):end-round(0.07*norm_size(1)), :);
    writeVideo(VID, fig_frame);
end
close (VID)

disp('ALL DONE')
end