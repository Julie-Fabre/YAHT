function bezier_gui()

% Create the main figure
hFig = figure('Position', [100, 100, 800, 600], 'Name', '3D Image Bezier GUI', 'NumberTitle', 'Off', 'MenuBar', 'None', 'Toolbar', 'None');

% Sample 3D data (Replace this with your image stack)
imgStack = rand(100, 100, 10); % 10 slices of 100x100 images

% Parameters
% Parameters
currentSlice = 1;
sliceControlPoints = cell(size(imgStack, 3), 1); % Store control points for each slice
controlPoints = sliceControlPoints{currentSlice}; % This is the active set for the current slice
hHighlighted = [];
pointSelected = false;


% GUI Controls
uicontrol('Style', 'text', 'Position', [10, 560, 100, 20], 'String', 'Slice:');
uicontrol('Style', 'slider', 'Position', [60, 560, 300, 20], 'Min', 1, 'Max', size(imgStack, 3), 'Value', ...
    currentSlice, 'SliderStep', [1 / (size(imgStack, 3) - 1), 1 / (size(imgStack, 3) - 1)], 'Callback', @(src, event)updateSlice(src));

hAxes = axes('Parent', hFig, 'Position', [0.1, 0.1, 0.7, 0.7]);
imagesc(imgStack(:, :, currentSlice), 'Parent', hAxes);
axis image off;
colormap gray;

% Initialize plots for control points and Bezier curve
hold on;
hPointsPlot = plot(hAxes, NaN, NaN, 'ro', 'ButtonDownFcn', @startDragFcn, 'Tag', 'draggable'); % initialize an empty plot for control points
hCurvePlot = plot(hAxes, NaN, NaN, 'r', 'LineWidth', 2); % initialize an empty plot for Bezier curve
hHighlighted = scatter(hAxes, NaN, NaN, 'b*', 'ButtonDownFcn', @startDragFcn, 'Tag', 'draggable');
hold off;
hListbox = uicontrol('Style', 'listbox', 'Position', [650, 250, 130, 300], 'String', {}, 'Min', 0, 'Max', 2, 'Callback', @(src, event)highlightSelectedPoint());
uicontrol('Style', 'pushbutton', 'Position', [650, 210, 130, 30], 'String', 'Delete Point', 'Callback', @(src, event)deleteSelectedPoint());


set(hFig, 'WindowButtonDownFcn', @(src, event)selectPoint(src));

    function updateSlice(src)
        currentSlice = round(get(src, 'Value'));
        controlPoints = sliceControlPoints{currentSlice}; % Load the control points for the current slice
        if isempty(controlPoints)
            controlPoints = [];
        end
        refreshDisplay();
    end

    function selectPoint(~)
        if ~pointSelected
            pt = get(hAxes, 'CurrentPoint');
            controlPoints(end+1, :) = pt(1, 1:2);

            sliceControlPoints{currentSlice} = controlPoints; % Store the updated control points
            updatePlots();
            updateListbox();
        end
    end

    function updateListbox()
        pointStrings = cell(size(controlPoints, 1), 1);
        for i = 1:size(controlPoints, 1)
            pointStrings{i} = ['Point ', num2str(i), ': (', num2str(controlPoints(i, 1)), ', ', num2str(controlPoints(i, 2)), ')'];
        end
        set(hListbox, 'String', pointStrings);
    end


    function deleteSelectedPoint()
        selectedIdx = get(hListbox, 'Value');
        if ~isempty(selectedIdx)
            controlPoints(selectedIdx, :) = [];
            sliceControlPoints{currentSlice} = controlPoints; % Store the updated control points
            removeHighlighted(); % Remove highlighted point
            updatePlots();
            updateListbox();
        end
    end

    function highlightSelectedPoint()
        removeHighlighted();
        selectedIdx = get(hListbox, 'Value');
        if ~isempty(selectedIdx)
            set(hHighlighted, 'XData', controlPoints(selectedIdx, 1), 'YData', controlPoints(selectedIdx, 2));
            pointSelected = true; % Indicate that a point is currently selected
        end
    end

    function removeHighlighted()
        if ishandle(hHighlighted)
            set(hHighlighted, 'XData', NaN, 'YData', NaN);
        end
    end


    function startDragFcn(~, ~)
        selectedIdx = get(hListbox, 'Value');
        if ~isempty(selectedIdx)
            set(hFig, 'WindowButtonMotionFcn', @draggingFcn);
            set(hFig, 'WindowButtonUpFcn', @stopDragFcn);
        end
    end

    function draggingFcn(~,~)
        pt = get(hAxes, 'CurrentPoint');
        selectedIdx = get(hListbox, 'Value');
        controlPoints(selectedIdx, :) = pt(1, 1:2);
        
        sliceControlPoints{currentSlice} = controlPoints; % Store the updated control points

        set(hHighlighted, 'XData', pt(1, 1), 'YData', pt(1, 2));
        updatePlots();
        updateListbox();
    end

    function stopDragFcn(~, ~)
        set(hFig, 'WindowButtonMotionFcn', '');
        set(hFig, 'WindowButtonUpFcn', '');
        highlightSelectedPoint(); % Re-highlight to reflect the position change
        pointSelected = false; % Reset the flag since we're no longer selecting a point
    end

    function refreshDisplay()
        imagesc(imgStack(:, :, currentSlice), 'Parent', hAxes);
        axis image off;
        colormap gray;
        title(['Slice: ', num2str(currentSlice)]);
        updatePlots();
    end

    function updatePlots()
        % Update the XData and YData of the control points and Bezier curve plots
        if ~isempty(controlPoints)
            set(hPointsPlot, 'XData', controlPoints(:, 1), 'YData', controlPoints(:, 2));

            if size(controlPoints, 1) > 1
                t = linspace(0, 1, 100);
                B = bezier(controlPoints', t);
                set(hCurvePlot, 'XData', B(1, :), 'YData', B(2, :));
            else
                set(hCurvePlot, 'XData', NaN, 'YData', NaN);
            end
        else
            set(hPointsPlot, 'XData', NaN, 'YData', NaN);
            set(hCurvePlot, 'XData', NaN, 'YData', NaN);
        end
    end

    function B = bezier(controlPoints, t)
        n = size(controlPoints, 2) - 1;
        B = zeros(2, length(t));
        for i = 0:n
            B = B + nchoosek(n, i) * (1 - t).^(n - i) .* t.^i .* controlPoints(:, i+1);
        end
    end

end
