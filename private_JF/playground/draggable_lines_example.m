function draggable_lines_example
    fig = figure;
    ax = axes('Parent', fig);
    hold(ax, 'on');

    % Example data for boundaries and plot
    trajectory_area_boundaries = [1, 2, 3];
    probe_areas_ax = ax;

    % Draw boundary lines and make them draggable
    boundary_lines = gobjects(length(trajectory_area_boundaries), 1);
    for curr_boundary = 1:length(trajectory_area_boundaries)
        boundary_lines(curr_boundary) = line(probe_areas_ax, [-13.5, 1.5], ...
            repmat((trajectory_area_boundaries(curr_boundary) - 0.5) * 25, 1, 2), 'color', 'b', 'linewidth', 2);
        draggable_line(boundary_lines(curr_boundary), ax, fig);
    end

    function draggable_line(hLine, ax, fig)
        set(hLine, 'ButtonDownFcn', @startDragFcn);
        
        function startDragFcn(src, ~)
            src.UserData = get(ax, 'CurrentPoint');
            set(fig, 'WindowButtonMotionFcn', {@draggingFcn, src});
            set(fig, 'WindowButtonUpFcn', @stopDragFcn);
        end
        
        function draggingFcn(~, ~, src)
            curr_pt = get(ax, 'CurrentPoint');
            orig_pt = src.UserData;
            
            dy = curr_pt(1, 2) - orig_pt(1, 2);
            ydata = get(src, 'YData');
            ydata = ydata + dy;
            set(src, 'YData', ydata);
            
            src.UserData = curr_pt;
        end
        
        function stopDragFcn(~, ~)
            set(fig, 'WindowButtonMotionFcn', '');
            set(fig, 'WindowButtonUpFcn', '');
        end
    end


end
