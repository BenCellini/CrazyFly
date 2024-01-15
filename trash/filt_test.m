close all ; clc

vid = squeeze(vidData);
dim = size(vid);

n = 5;
Fc = 25;
Fs = 200;
filtvid = filtfilt_vid(vid, n, Fc, Fs);

%%
close all ; clc

frame = vid(:,:,1);

M = cell(1,1);
M{1} = frame;
M{end+1} = imsharpen(M{end});
M{end+1} = 1.5*medfilt2(M{end},[9 9]);
M{end+1} = imsharpen(M{end});
M{end+1} = imgaussfilt(M{end});
M{end+1} = 1.5*M{end};
M{end+1} = imsharpen(M{end});
M{end+1} = medfilt2(M{end},[9 9]);
% M{end+1} = imbinarize(M{end}, 'adaptive', 'Sensitivity', 'ForegroundPolarity', 'bright');
montage(M)

%%
close all ; clc

tic
thresh = 35;
win = 5;
step = 3;
c_step = ceil(dim(2) / win);
r_step = ceil(dim(1) / win);
test_frame = M{end};
imshow(test_frame) ; hold on
mask_init = false(size(test_frame));
roi_block = nan(r_step, c_step);
rb = 1;
new_frame = test_frame;
for r = 1:step:dim(1)
    cb = 1;
    for c = 1:step:dim(2)
        x = [c, c + win];
     	y = [r, r + win];
        
        x(x > dim(2)) = dim(2);
        y(y > dim(1)) = dim(1);
        
        mask = mask_init;
        x_span = x(1):x(2);
        y_span = y(1):y(2);
        mask(y_span,x_span) = true;
        
        roi = test_frame(mask);
        roi_int = mean(roi);
        roi_block(rb,cb) = roi_int;
        
        %new_frame(y_span,x_span) = roi_int;
        if roi_int < thresh
            new_frame(y_span,x_span) = 0;
        else
            %new_frame(y_span,x_span) = 0;
        end
                        
        %imshow(mask)
        %pause(0.1)
        
%         h = rectangle('Position', [c r win win], 'EdgeColor', 'c', 'FaceColor', 'r');
%         pause(0.1)
%         delete(h)

        cb = cb + 1;
    end
    rb = rb + 1;
end
toc
imshow(new_frame)


% filt_frame = filtvid(:,:,1);
% 
% 
% M{2} = 2*filt_frame;

