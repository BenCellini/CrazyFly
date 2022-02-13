function [angle, h] = set_angle(frame, color)
%% set_angle: manually set the angle of an object in a frame
%
%   INPUT:
%       frame 	: 2D image matrix
%       color 	: color of line
%
%   OUTPUT:
%       angle	: manually defined angle
%       h       : line graphic object
%

if nargin < 2
    color = [1 0 0];
end

clf ; hold on
imshow(frame) ; title('Press escape when done')
roi = drawline(gca, 'Color', color);
addlistener(roi,'ROIMoved',@allevents);
wait(roi)

global h
angle = calculate_angle(h);
close

end

function allevents(src,evt)
    global h
    evname = evt.EventName;
    switch(evname)
        case{'ROIMoved'}
            h = src.Position;
            angle = calculate_angle(h);
            disp(['Angle: ' num2str(angle)]);
    end
end

function angle = calculate_angle(pos)
    dx = pos(1,1) - pos(2,1);
    dy = pos(1,2) - pos(2,2);
    angle = rad2deg(atan2(dy, dx)) + 90;
end