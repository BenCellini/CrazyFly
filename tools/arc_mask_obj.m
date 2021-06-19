classdef arc_mask_obj
    % arc_mask_obj: make a UI arc mask
    %   
    
    properties
        H       % graphics handles structure
        sym     % symmetric boolean
    end
    
    properties (SetAccess = public)
        iter            % iteration of mask
        type            % global or local
        move            % movable points and lines
        move_points     % position of movable points
        inner           % inner
        outer           % outer
        global_ang      % global angle [°]
        radius          % radius structure
        span            % angular span of mask [°]
        points          % mask movable points
        left_angle      % left angle in global [°]
        right_angle     % right angle in global [°]
        angles          % angles between left & right
        init            % intial angle settings
        axis
        axis_reverse
        init_line
        area_points
        ROI
        image
        image_sz
        frame
        frame_bw
        bw_thresh       
        y
        x
        lw              % line width
        color           % color
        mrk_sz          % marker size
        patch_h         % patch
    end
    
    methods
        function mask = arc_mask_obj(frame, rot, global_ang, R, span, mask_color, iter, type)
            % arc_mask_obj: Construct an instance of this class
            %   
                        
            fig = gcf;
            
            if nargin < 8
                type = 1;
                if nargin < 7
                    iter = 1;
                end
            end
            mask.type = type;
            mask.iter = iter;
            
            switch mask.type
                case 1
                    mask.init.color = [1 0 0];
                    axis_scale = 2;
                case 2
                    mask.init.color = 'k';
                    axis_scale = 0.2;
                otherwise
                    error('"type" must be 1 or 2')
            end
            
            % Defaults
            mask.lw = 0.1;
            mask.color = colorspec2rgb(mask_color); % mask color
            mask.sym = false;
            
            % Get frame to set mask
            mask.image = frame;
            mask.image_sz = size(mask.image);

            % Intiial mask settings
            mask.global_ang = global_ang; % body angle (mask reference frame) [°]
            mask.radius.inner = R(1); % inner mask radius
            mask.radius.outer = R(2); % outer mask radius
            mask.radius.axis = 0.5*mask.radius.outer; % axis radius
            mask.span = span; % angular span of mask [+,-] [°]

            % Define rotation point
            hold on
            mask.move.rot = drawpoint('Position', rot, 'Color', mask.color);
            mask.mrk_sz = 0.7 * mask.move.rot.MarkerSize;
            mask.move.rot.MarkerSize = mask.mrk_sz;

            % Define mask angles in global frame
            mask.left_angle = mask.global_ang - mask.span(1); % left angle in global [°]
            mask.right_angle = mask.global_ang + mask.span(2); % right angle in global [°]
            mask.angles = (mask.left_angle:1:mask.right_angle)'; % angles in between edges of mask

            % Define initial reference angle & radius
            mask.init.radius = mean([mask.radius.inner, mask.radius.outer]); % radius of initial refercne angle
            mask.init.angle = 0; % initial reference angle [°]
            
            % Define arc coordinates
            inner_cent_pos = mask.move.rot.Position + ...
                mask.radius.inner*[sind(mask.global_ang) , -cosd(mask.global_ang)];
            outer_cent_pos = mask.move.rot.Position + ...
                mask.radius.outer*[sind(mask.global_ang) , -cosd(mask.global_ang)];
            outer_left_pos = mask.move.rot.Position + ...
                mask.radius.outer*[sind(mask.left_angle) , -cosd(mask.left_angle)];
            outer_right_pos = mask.move.rot.Position + ...
                mask.radius.outer*[sind(mask.right_angle) , -cosd(mask.right_angle)];
            axis_pos = mask.move.rot.Position - ...
                axis_scale*mask.radius.axis*[sind(mask.global_ang) , -cosd(mask.global_ang)];

            % Point to define initial reference angle
            init_pos = mask.move.rot.Position + ...
                mask.init.radius*[sind(mask.global_ang + mask.init.angle) , ...
                -cosd(mask.global_ang + mask.init.angle)];

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
            mask.move.init = drawpoint('Position', init_pos, ...
                'Color', mask.init.color, 'MarkerSize', ceil(1.2*mask.mrk_sz));

            % UI Controls
            if mask.type == 1
                fig_h = fig.Position(4);
                mask.H.sym = uicontrol('style','checkbox','units','pixels', 'Value', mask.sym, ...
                                'position',[10,0.055*fig_h,70,15],'string','symmetric');
                set(mask.H.sym, 'callback', @(src, event)sym_call(mask ,src, event));

                mask.H.global_ang = uicontrol('style','edit','units','pixels', ...
                                'position',[10,0.15*fig_h,50,15],'string','0');
                set(mask.H.global_ang, 'callback', @(src, event) global_call(mask, src, event), ...
                    'BackgroundColor', [0 1 0]);

                mask.H.init = uicontrol('style','edit','units','pixels', ...
                                'position',[10,0.1*fig_h,50,15],'string','0');
                set(mask.H.init, 'callback', @(src, event) init_call(mask, src, event), ...
                    'BackgroundColor', mask.init.color);

                mask.H.done = uicontrol('style','pushbutton','units','pixels', 'Value', false, ...
                                'position',[10,0.01*fig_h,50,15],'string','Done');
                set(mask.H.done, 'callback', @(src, event) done_call(mask, src, event));
            end

            % Draw mask
            mask.patch_h = [];
            mask.axis = [];
            mask.axis_reverse = [];
            mask.init_line = [];
            mask = draw_mask(mask);
            mask = get_ROI(mask);
            
            % Store mask in GUI
            gdata = guidata(gcf);
            gdata{mask.iter} = mask;
            guidata(fig, gdata)

            % Listeners to move mask
            addlistener(mask.move.outer_C, 'MovingROI', @(src,evnt)outer_C(mask,src,evnt));
            addlistener(mask.move.outer_C, 'ROIMoved', @(src,evnt)outer_C(mask,src,evnt));
            addlistener(mask.move.outer_L, 'MovingROI', @(src,evnt)outer_L(mask,src,evnt));
            addlistener(mask.move.outer_L, 'ROIMoved', @(src,evnt)outer_L(mask,src,evnt));
            addlistener(mask.move.outer_R, 'MovingROI', @(src,evnt)outer_R(mask,src,evnt));
            addlistener(mask.move.outer_R, 'ROIMoved', @(src,evnt)outer_R(mask,src,evnt));
            addlistener(mask.move.inner_C, 'MovingROI', @(src,evnt)inner_C(mask,src,evnt));
            addlistener(mask.move.inner_C, 'ROIMoved', @(src,evnt)inner_C(mask,src,evnt));
            addlistener(mask.move.rot, 'MovingROI', @(src,evnt)rot_point(mask,src,evnt));
            addlistener(mask.move.rot, 'ROIMoved', @(src,evnt)rot_point(mask,src,evnt));
            addlistener(mask.move.axis, 'MovingROI', @(src,evnt)axis_h(mask,src,evnt));
            addlistener(mask.move.axis, 'ROIMoved', @(src,evnt)axis_h(mask,src,evnt));
            addlistener(mask.move.init, 'MovingROI', @(src,evnt)init_h(mask,src,evnt));
            addlistener(mask.move.init, 'ROIMoved', @(src,evnt)init_h(mask,src,evnt));
        end
        
        function mask = draw_mask(mask)
            % Define new mask angles
            %mask.span(mask.span < 0) = mask.span(mask.span < 0) + 360;
            mask.angles = (mask.global_ang + ( -mask.span(1):mask.span(2) ))';
            mask.left_angle = mask.angles(1);
            mask.right_angle = mask.angles(end);

            % Create the new arc mask point map
            mask.inner.pts = mask.radius.inner*[sind(mask.angles) , -cosd(mask.angles)];
            mask.outer.pts = mask.radius.outer*[sind(mask.angles) , -cosd(mask.angles)];
            mask.points = mask.move.rot.Position + ...
                [mask.inner.pts ; flipud(mask.outer.pts) ; mask.inner.pts(1,:)];

            % Constrain points
            mask.move.inner_C.Position = mask.move.rot.Position + ...
                mask.radius.inner*[sind(mask.global_ang) , -cosd(mask.global_ang)];
            mask.move.outer_C.Position = mask.move.rot.Position + ...
                mask.radius.outer*[sind(mask.global_ang) , -cosd(mask.global_ang)];
            mask.move.outer_L.Position = mask.move.rot.Position + ...
                mask.radius.outer*[sind(mask.left_angle) , -cosd(mask.left_angle)];
            mask.move.outer_R.Position = mask.move.rot.Position + ...
                mask.radius.outer*[sind(mask.right_angle) , -cosd(mask.right_angle)];
            mask.move.axis.Position = mask.move.rot.Position - ...
                mask.radius.axis*[sind(mask.global_ang) , -cosd(mask.global_ang)];
            mask.move.init.Position = mask.move.rot.Position + ...
                mask.init.radius*[sind(mask.global_ang + mask.init.angle) , ...
                -cosd(mask.global_ang + mask.init.angle)];

            % Annotate the mask ROI reigon
            delete(mask.patch_h)
            delete([mask.axis mask.axis_reverse mask.init_line])
            mask.patch_h = patch(mask.points(:,1), mask.points(:,2), ...
                mask.color, 'FaceAlpha', 0.1, 'EdgeColor', mask.color, 'LineWidth', 0.75);
            mask.axis = plot([mask.move.rot.Position(1) mask.move.outer_C.Position(1)], ...
                [mask.move.rot.Position(2) mask.move.outer_C.Position(2)], ...
                 '--', 'Color', mask.color, 'LineWidth', 0.25);
            mask.axis_reverse = plot([mask.move.axis.Position(1) mask.move.rot.Position(1)], ...
                [mask.move.axis.Position(2) mask.move.rot.Position(2)], ...
                 '--', 'Color', mask.color, 'LineWidth', 0.25);
            mask.init_line = plot([mask.move.rot.Position(1) mask.move.init.Position(1)], ...
                [mask.move.rot.Position(2) mask.move.init.Position(2)], ...
                 '-', 'Color', mask.init.color, 'LineWidth', 0.25);

            % Store the points
            mask.move_points = structfun(@(x) x.Position, mask.move, 'UniformOutput', false);

            % Update global_ang heading display
            mask.H.global_ang.String = string(mask.global_ang);
            mask.H.init.String = string(mask.init.angle);
            
            % Store mask in GUI
            gdata = guidata(gcf);
            gdata{mask.iter} = mask;
            guidata(gcf, gdata)
        end

        function mask = get_ROI(mask, frame, showplot)
            if nargin < 3
               showplot= false; 
               if nargin < 2
                  frame = mask.image; 
               end
            end
            
            if isempty(frame)
                frame = mask.image; 
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
            elseif nargin == 2
                mask.frame = frame;
            end        
            mask.frame = frame;
            mask.bw_thresh = graythresh(mask.frame(mask.area_points));
            mask.frame(~mask.area_points) = 0;
            mask.y = round(mask.ROI(3):mask.ROI(4));
            mask.x = round(mask.ROI(1):mask.ROI(2));
            
            mask.frame = mask.frame(mask.y, mask.x);
            mask.frame_bw = imbinarize(mask.frame, mask.bw_thresh);

            if showplot
                fig = figure;
                set(fig, 'Name','ROI')
                imshow(mask.frame_bw)
            end
            
          	% Store mask in GUI
            gdata = guidata(gcf);
            gdata{mask.iter} = mask;
            guidata(gcf, gdata)
        end
        
        function rot_point(obj,src,evt)
            mask = guidata(src);
            mask = mask{obj.iter};
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

        function axis_h(obj,src,evt)
            mask = guidata(src);
            mask = mask{obj.iter};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = evt.CurrentPosition - mask.move.rot.Position;
                    mask.radius.axis = sqrt(shift(1)^(2) + shift(2)^(2));
                    mask.global_ang = rad2deg(atan2(shift(2), shift(1))) - 90;
                    mask = draw_mask(mask);
                case{'ROIMoved'}
                    mask = get_ROI(mask);
            end
        end
        
        function outer_C(obj,src,evt)
            mask = guidata(src);
            mask = mask{obj.iter};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = evt.CurrentPosition - mask.move.rot.Position;
                    mask.radius.outer = sqrt(shift(1)^(2) + shift(2)^(2));
                    %mask.global_ang = rad2deg(atan2(shift(2), shift(1))) + 90;
                    mask = draw_mask(mask);
                case{'ROIMoved'}
                    mask = get_ROI(mask);
            end
        end
                
        function outer_L(obj,src,evt)
            mask = guidata(src);
            mask = mask{obj.iter};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = mask.move.rot.Position - evt.CurrentPosition;
                    mask.span(1) = rad2deg(atan2(shift(1), shift(2))) + mask.global_ang;
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

        function outer_R(obj,src,evt)
            mask = guidata(src);
            mask = mask{obj.iter};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = mask.move.rot.Position - evt.CurrentPosition;
                    mask.span(2) = rad2deg(atan2(-shift(1), shift(2))) - mask.global_ang;
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
        
        function inner_C(obj,src,evt)
            mask = guidata(src);
            mask = mask{obj.iter};
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
        
        function init_h(obj,src,evt)
            mask = guidata(src);
            mask = mask{obj.iter};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = mask.move.rot.Position - evt.CurrentPosition;
                    mask.init.angle = rad2deg(atan2(-shift(1), shift(2))) - mask.global_ang;
                    mask.init.radius = sqrt(shift(1)^(2) + shift(2)^(2));
                    mask = draw_mask(mask);
                case{'ROIMoved'}
                    mask = get_ROI(mask);
            end
        end
        
        function sym_call(obj,src,event)
            mask = guidata(src);
            mask = mask{obj.iter};
            mask.sym = mask.H.sym.Value;
            
            gdata = guidata(gcf);
            gdata{mask.iter} = mask;
            guidata(gcf, gdata)
        end
        
        function global_call(obj,src,event)
            mask = guidata(src);
            mask = mask{obj.iter};
            
            mask.global_ang = str2double(mask.H.global_ang.String);
            mask = draw_mask(mask);
        end
        
        function init_call(obj,src,event)
            mask = guidata(src);
            mask = mask{obj.iter};
            
            mask.init.angle = str2double(mask.H.init.String);
            mask = draw_mask(mask);
        end
        
        function done_call(obj,src,event)
            close(gcf)
        end 
    end
    
end

