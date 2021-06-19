classdef fly_mask
    % fly_mask: make a UI arc mask
    %   
    
    properties
        H       % graphics handles structure
        sym     % symmetric boolean
    end
    
    properties (SetAccess = public)
        iter            % iteration of mask
        move            % movable points and lines
        lwing
        rwing
        guiI
        axis
        cent
        radius
        global_ang
        wing_cent
    end
    
    methods
        function mask = fly_mask(frame, neck, global_ang)
            % arc_mask_obj: Construct an instance of this class
            %   
            
            mask.guiI = 3; % GUi data index
            mask.sym = true; % make mask symmetric
            mask.global_ang = global_ang; % global mask angle
            
            mask.radius.axis = 55; % body axis
            mask.radius.xoffset = 20; % wing joint distance from neck (x)
            mask.radius.yoffset = 35; % wing joint distance from neck (y)
            
            abdomen = neck + mask.radius.axis* ... % abdomen point
            [-sind(mask.global_ang ) , cosd(mask.global_ang )];
            
            mask.wing_cent = neck + mask.radius.yoffset*... % wing center point
                [-sind(mask.global_ang) , cosd(mask.global_ang)];
                        
            lwing_rot = mask.wing_cent + mask.radius.xoffset* ... % left wing joint point
                [sind(mask.global_ang - 90 ) , -cosd(mask.global_ang - 90)];
            rwing_rot = mask.wing_cent + mask.radius.xoffset* ...  % right wing joint point
                [sind(mask.global_ang + 90 ) , -cosd(mask.global_ang + 90 )];
            
            % Make the movable mask ROI points
            mask.move.neck = drawpoint('Position', neck, ...
                'Color', 'c', 'MarkerSize', 5);
            mask.move.abdomen = drawpoint('Position', abdomen, ...
                'Color', 'm', 'MarkerSize', 5);
                          
            mask.axis = [];
            mask.cent = [];
            
            % Make arc masks for each wing
            mask.lwing = arc_mask_obj(frame, lwing_rot, ...
                mask.global_ang-90, [50 100], [20 50], 'g', 1, 2);
            mask.rwing = arc_mask_obj(frame, rwing_rot, ...
                mask.global_ang+90, [50 100], [50 20], 'b', 2, 2);
            
            mask.lwing.init.color = 'none';
            mask.rwing.init.color = 'none';
            
            % Draw the mask
            mask = draw(mask);
            
            % Store mask in GUI
            gdata = guidata(gcf);
            gdata{mask.guiI} = mask;
            guidata(gcf, gdata)
            
            % Listeners to move mask
            addlistener(mask.move.neck, 'MovingROI', @(src,evnt)neck_h(mask,src,evnt));
            addlistener(mask.move.neck, 'ROIMoved', @(src,evnt)neck_h(mask,src,evnt));
            addlistener(mask.move.abdomen, 'MovingROI', @(src,evnt)abdomen_h(mask,src,evnt));
            addlistener(mask.move.abdomen, 'ROIMoved', @(src,evnt)abdomen_h(mask,src,evnt));  
            addlistener(mask.lwing.move.rot, 'MovingROI', @(src,evnt)lwing_h(mask,src,evnt));
            addlistener(mask.lwing.move.rot, 'ROIMoved', @(src,evnt)lwing_h(mask,src,evnt));
            addlistener(mask.rwing.move.rot, 'MovingROI', @(src,evnt)rwing_h(mask,src,evnt));
            addlistener(mask.rwing.move.rot, 'ROIMoved', @(src,evnt)rwing_h(mask,src,evnt));
        end
        
        function mask = draw(mask)
            mask.lwing.global_ang = mask.global_ang - 90;
            mask.rwing.global_ang = mask.global_ang + 90;
            mask.wing_cent = mask.move.neck.Position + mask.radius.yoffset*...
                [-sind(mask.global_ang ) , cosd(mask.global_ang )];

            mask.lwing.move.rot.Position = mask.wing_cent + mask.radius.xoffset* ....
                [sind(mask.global_ang -90) , -cosd(mask.global_ang - 90)];
            
            mask.rwing.move.rot.Position = mask.wing_cent + mask.radius.xoffset* ....
                [sind(mask.global_ang + 90) , -cosd(mask.global_ang + 90)];
            
            draw_mask(mask.lwing);
            draw_mask(mask.rwing);
            delete([mask.axis mask.cent])
            mask.axis = plot([mask.move.neck.Position(1) mask.move.abdomen.Position(1)], ...
                [mask.move.neck.Position(2) mask.move.abdomen.Position(2)], ...
                 '--', 'Color', 'c', 'LineWidth', 0.25);
            mask.cent = plot(mask.wing_cent(1), mask.wing_cent(2), ...
                 '*', 'Color', 'y', 'MarkerSize', 3);
             
            % Store mask in GUI
            gdata = guidata(gcf);
            gdata{mask.guiI} = mask;
            guidata(gcf, gdata)
        end

        function mask = get_mask(mask, frame)
            if nargin < 2
                frame = [];
            end

            mask.lwing = get_ROI(mask.lwing, frame);
            mask.rwing = get_ROI(mask.rwing, frame);
            
            global mask

%             % Store mask in GUI
%             gdata = guidata(gcf);
%             gdata{mask.guiI} = mask;
%             guidata(gcf, gdata)
        end
        
        function neck_h(obj,src,evt)
            temp = guidata(src);
            mask = temp{obj.guiI};
            mask.lwing = temp{1};
            mask.rwing = temp{2};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = evt.CurrentPosition - evt.PreviousPosition;
                    mask.lwing.move.rot.Position = mask.lwing.move.rot.Position + shift;
                    mask.rwing.move.rot.Position = mask.rwing.move.rot.Position + shift;
                    mask.move.abdomen.Position = mask.move.abdomen.Position + shift;
                    draw(mask);
                case{'ROIMoved'}
                    get_mask(mask);
            end
        end
        
        function abdomen_h(obj,src,evt)
            temp = guidata(src);
            mask = temp{obj.guiI};
            mask.lwing = temp{1};
            mask.rwing = temp{2};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = evt.CurrentPosition - mask.move.neck.Position;
                    mask.radius.axis = sqrt(shift(1)^(2) + shift(2)^(2));
                    mask.global_ang = rad2deg(atan2(shift(2), shift(1))) - 90;
                    draw(mask);
                case{'ROIMoved'}
                    get_mask(mask);
            end
        end
        
        function lwing_h(obj,src,evt)
            temp = guidata(src);
            mask = temp{obj.guiI};
            mask.lwing = temp{1};
            mask.rwing = temp{2};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = mask.move.neck.Position - evt.CurrentPosition;
                    r = sqrt(shift(1)^(2) + shift(2)^(2));
                    ang = rad2deg(atan2(shift(2), shift(1))) -  mask.global_ang;
                    mask.radius.xoffset = r*cosd(ang);
                    mask.radius.yoffset = -r*sind(ang);
                    draw(mask);
                case{'ROIMoved'}
                    get_mask(mask);
            end
        end
        
        function rwing_h(obj,src,evt)
            temp = guidata(src);
            mask = temp{obj.guiI};
            mask.lwing = temp{1};
            mask.rwing = temp{2};
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
                    shift = mask.move.neck.Position - evt.CurrentPosition;
                    r = sqrt(shift(1)^(2) + shift(2)^(2));
                    ang = rad2deg(atan2(shift(2), shift(1))) -  mask.global_ang;
                    mask.radius.xoffset = -r*cosd(ang);
                    mask.radius.yoffset = -r*sind(ang);
                    draw(mask);
                case{'ROIMoved'}
                    get_mask(mask);
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

