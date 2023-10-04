function bd_drawProbes_bezierCurves(tv, av, st, registeredImage, outputDir, screenToUse)
% based on:
% AP_get_probe_histology(tv,av,st,slice_im_path)
%
% QQ: fine-tune button locations: they are hard-coded and only work work for 1 screen type

%% create the GUI figure
% Initialize guidata
gui_data = struct;
gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;

% Query number of probes from user
gui_data.n_probes = str2num(cell2mat(inputdlg('How many probes?')));

% Load in slice images

gui_data.slice_im_path = outputDir;

for curr_slice = 1:size(registeredImage, 3)
    gui_data.slice_im{curr_slice} = registeredImage(:, :, curr_slice);
end

% Load corresponding CCF slices
ccf_slice_fn = [outputDir, filesep, '/manual/histology_ccf.mat'];
load(ccf_slice_fn);
gui_data.histology_ccf = histology_ccf;

% Load histology/CCF alignment
ccf_alignment_fn = [outputDir, filesep, '/manual/atlas2histology_tform.mat'];
load(ccf_alignment_fn);
gui_data.histology_ccf_alignment = atlas2histology_tform;

% Warp area labels by histology alignment
gui_data.histology_aligned_av_slices = cell(length(gui_data.slice_im), 1);
for curr_slice = 1:length(gui_data.histology_ccf)
    curr_av_slice = squeeze(gui_data.histology_ccf(curr_slice).av_slices);
    curr_av_slice(isnan(curr_av_slice)) = 1;
    curr_slice_im = gui_data.slice_im{curr_slice};

    tform = affine2d;
    tform.T = gui_data.histology_ccf_alignment{curr_slice};
    tform_size = imref2d([size(curr_slice_im, 1), size(curr_slice_im, 2)]);
    gui_data.histology_aligned_av_slices{curr_slice} = ...
        imwarp(curr_av_slice, tform, 'nearest', 'OutputView', tform_size);
end

% Create figure, set button functions
gui_fig = figure('KeyPressFcn', @keypress, 'Color', 'k');

SCRSZ = screensize(screenToUse); %Get user's screen size
gui_data.SCRSZ = SCRSZ;
set(gui_fig, 'Position', SCRSZ);
screenPortrait = SCRSZ(4) > SCRSZ(3);
if screenPortrait
    gui_data.histology_ax = axes('Position', [0.05, 0.5, 0.9, 0.45], 'YDir', 'reverse');
    gui_button_position1 = [SCRSZ(3) * 0.1, SCRSZ(4) * 0.45, 100, 40];
    gui_button_position2 = [SCRSZ(3) * 0.1, SCRSZ(4) * 0.9, 100, 40];
else
    gui_data.histology_ax = axes('Position', [0.1, 0.1, 0.5, 0.5], 'YDir', 'reverse');
    gui_button_position1 = [SCRSZ(3) * 0.1, SCRSZ(4) * 0.45, 100, 40]; %QQ change
    gui_button_position2 = [SCRSZ(3) * 0.1, SCRSZ(4) * 0.45, 100, 40]; %QQ change

end
gui_data.gui_button_position = gui_button_position1;

% auto contrast/brightness button
gui_data.auto_contrast_btn = uicontrol('Style', 'pushbutton', ...
    'String', '<HTML><center><FONT color="white"><b>Auto brightness/ contrast</b></Font>', ...
    'Position', gui_button_position2+[0, 0, 0, 20], ...
    'BackgroundColor', rgb('DeepPink'), ...
    'CallBack', @(varargin) autoContrastButtonPushed(gui_fig));

% brightness slider
gui_data.brightness_beta = 0;
gui_data.brightness_slider = uicontrol('Style', 'slider', ...
    'String', 'Brightness', ...
    'Position', gui_button_position2+[150, 0, 50, -20], ...
    'BackgroundColor', rgb('Black'), ...
    'Min', -100, 'Max', 100, 'Value', gui_data.brightness_beta, ...
    'CallBack', @(varargin) brightnessButtonPushed(gui_fig));

gui_data.brightness_text = uicontrol('Style', 'text', ...
    'String', 'Brightness %', ...
    'Position', gui_button_position2+[150, 20, 50, -20], ...
    'BackgroundColor', rgb('Black'), ...
    'ForegroundColor', rgb('White'));

% contrast slider
gui_data.contrast_alpha = 1;
gui_data.contrast_slider = uicontrol('Style', 'slider', ...
    'String', 'Contrast', ...
    'Position', gui_button_position2+[400, 0, 50, -20], ...
    'BackgroundColor', rgb('Black'), ...
    'Min', 0, 'Max', 1, 'Value', gui_data.contrast_alpha, ...
    'CallBack', @(varargin) contrastButtonPushed(gui_fig));

gui_data.contrast_text = uicontrol('Style', 'text', ...
    'String', 'Contrast %', ...
    'Position', gui_button_position2+[400, 20, 50, -20], ...
    'BackgroundColor', rgb('Black'), ...
    'ForegroundColor', rgb('White'));

% point size slider
gui_data.point_size = 20;
gui_data.point_size_slider = uicontrol('Style', 'slider', ...
    'String', 'Point size', ...
    'Position', gui_button_position2+[650, 0, 50, -20], ...
    'BackgroundColor', rgb('Black'), ...
    'Min', 1, 'Max', 100, 'Value', gui_data.point_size, ...
    'CallBack', @(varargin) pointSizeButtonPushed(gui_fig));

gui_data.point_size_text = uicontrol('Style', 'text', ...
    'String', 'Point size (px)', ...
    'Position', gui_button_position2+[650, 20, 50, -20], ...
    'BackgroundColor', rgb('Black'), ...
    'ForegroundColor', rgb('White'));

% save current datapoints
gui_data.save_btn = uicontrol('Style', 'pushbutton', ...
    'String', '<HTML><center><FONT color="white"><b>Save</b></Font>', ...
    'Position', [400, 200, 100, 50], ...
    'BackgroundColor', rgb('DarkOrange'), ...
    'CallBack', @(varargin) saveButtonPushed(gui_fig));

% load previous datapoints
gui_data.load_btn = uicontrol('Style', 'pushbutton', ...
    'String', '<HTML><center><FONT color="white"><b>Load</b></Font>', ...
    'Position', [400, 300, 100, 50], ...
    'BackgroundColor', rgb('DarkMagenta'), ...
    'CallBack', @(varargin) loadButtonPushed(gui_fig));


gui_data.curr_slice = 150;
gui_data.visibility = 1;
gui_data.fit_visibility = 0;

% Set up axis for histology image
hold on; colormap(gray); axis image off;
gui_data.slice_im{1}(gui_data.slice_im{1} > 1200) = 0;
img = imadjust(gui_data.slice_im{1}, [0.1, 0.8]);
img(img == 1) = 0;
gui_data.histology_im_h = image(img, ...
    'Parent', gui_data.histology_ax);
colormap(gray)

% Create title to write area in
gui_data.histology_ax_title = title(gui_data.histology_ax, '', 'FontSize', 14, 'Color', 'white');

% Initialize probe points -> old method
gui_data.probe_color = lines(gui_data.n_probes);
gui_data.probe_points_histology = cell(length(gui_data.slice_im), gui_data.n_probes);
for iProbe = 1:gui_data.n_probes
    gui_data.probe_points{iProbe} = plot(NaN, NaN, 'ro'); %, 'ButtonDownFcn', @(src, event)startDragFcn(src, event, gui_fig), 'Tag', 'draggable');
    gui_data.start_points_scatter{iProbe} = plot(NaN, NaN, 'rx');
    gui_data.stop_points_scatter{iProbe} = plot(NaN, NaN, 'rx');
    
end
gui_data.inflection_points_scatter = gobjects(gui_data.n_probes, 1);
gui_data.probe_fit_lines = gobjects(gui_data.n_probes, 1);
gui_data.stop_points = nan(gui_data.n_probes, 3);
gui_data.start_points = nan(gui_data.n_probes, 3);

% Initialize plots for control points and Bezier curve
points_table_position = [850, 120, 500, 300];
gui_data.bezier_control_points = cell(gui_data.n_probes, 1);
gui_data.bezier_curves = gobjects(gui_data.n_probes, 1);
gui_data.points_table = uicontrol('Parent', gui_fig, 'Style', 'listbox', 'Position', points_table_position, ...
    'String', {}, 'Callback', @(src, event)selectPoint(gui_fig, src, event));
gui_data.deletePointButton = uicontrol('Style', 'pushbutton', 'Position', ...
    [points_table_position(1,2) + [750, -50], 130, 30], 'String', 'Delete Point',...
    'BackgroundColor', rgb('Crimson'),...
    'Callback', @(src, event)deleteSelectedPoint(src, event, gui_fig));

% Start and stop points 
gui_data.probe_points_start = gobjects(gui_data.n_probes, 1);
gui_data.probe_points_start_position = nan(gui_data.n_probes, 3);
gui_data.probe_points_stop = gobjects(gui_data.n_probes, 1);
gui_data.probe_points_stop_position = nan(gui_data.n_probes, 3);
gui_data.startPt = uicontrol('Style', 'pushbutton', ...
        'String', 'Set starting point', ...
        'Position', [points_table_position(1,2)  + [880, -50], 130, 30], ...
        'BackgroundColor', rgb('PaleGreen'), ...
        'CallBack', @(varargin) startPtButtonPushed(gui_fig));

gui_data.stopPt = uicontrol('Style', 'pushbutton', ...
    'String', 'Set endpoint', ...
    'Position', [points_table_position(1,2) + [1070, -50], 130, 30], ...
    'BackgroundColor', rgb('LightCoral'), ...
    'CallBack', @(varargin) stopPtButtonPushed(gui_fig));


% Bezier Parameters
% gui_data.hHighlighted = [];
gui_data.pointSelected = false;
%gui_data.hHighlighted = scatter(NaN, NaN, 'b*', 'ButtonDownFcn', @(src, event)disp('Point clicked!'), 'Tag', 'draggable');
gui_data.hHighlighted = scatter(NaN, NaN, 40, 'b*', 'ButtonDownFcn', @(src, event)startDragFcn(src, event, gui_fig), 'Tag', 'draggable');
gui_data.selected_row = NaN;

% initialize probe buttons
nProbes_fit = floor((SCRSZ(4) - gui_button_position1(2))/50);
nCols = ceil(gui_data.n_probes/nProbes_fit);
colSpacing = (SCRSZ(3) - 100) / nCols;
gui_data.prev_spline_fit = repmat(0.1, gui_data.n_probes, 1);

% curr_probe
gui_data.curr_probe = 1;

% Create a cell array of probe names or identifiers
probes_list = arrayfun(@(x) ['Probe ', num2str(x)], 1:gui_data.n_probes, 'UniformOutput', false);
gui_data.probe_dropdown = uicontrol('Style', 'popupmenu', ...
    'String', probes_list, ...
    'Position', [points_table_position(1,2) + [730, 300], 500, 30], ...
    'Callback', @(src, event)updateProbe(gui_fig));

% Save gui_data
guidata(gui_fig, gui_data);

% Upload gui data
guidata(gui_fig, gui_data);

% Update the slice
update_slice(gui_fig);

end

%% key press
function keypress(gui_fig, eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key

    % left/right: move slice
    case 'leftarrow'
        gui_data.curr_slice = max(gui_data.curr_slice-1, 1);
        guidata(gui_fig, gui_data);
        update_slice(gui_fig);

    case 'rightarrow'
        gui_data.curr_slice = ...
            min(gui_data.curr_slice+1, length(gui_data.slice_im));
        guidata(gui_fig, gui_data);
        update_slice(gui_fig);


    case 'add' % 'add' is the numpad + key
        curr_probe = gui_data.probe_dropdown.Value;
        update_curr_probe(gui_fig, curr_probe)

    case '+'
        curr_probe = gui_data.probe_dropdown.Value;
        update_curr_probe(gui_fig, curr_probe)

    case 'insert' % add a probe
        addProbe(gui_fig)

    case 'uparrow'
        gui_data.curr_probe = min(gui_data.n_probes, gui_data.curr_probe+1);
        guidata(gui_fig, gui_data);
        updateProbe_slider(gui_fig, gui_data.curr_probe)

    case 'downarrow'
        gui_data.curr_probe = max(1, gui_data.curr_probe-1);
        guidata(gui_fig, gui_data);
        updateProbe_slider(gui_fig, gui_data.curr_probe)

    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?', 'Confirm exit', opts);
        if strcmp(user_confirm, 'Yes')

            saveButtonPushed(gui_fig)

            % Plot probe trajectories
            plot_probe(gui_data, probe_ccf);

        end
end

end

%% general functions: update probe, slice, create fit, probe buttons
function update_curr_probe(gui_fig, curr_probe)
% Get guidata
gui_data = guidata(gui_fig);
gui_data.curr_probe = curr_probe;

set(gui_data.histology_ax_title, 'String', ['Draw control points for probe ', num2str(curr_probe)]);
curr_point = drawpoint;

% Storing Bezier control points instead of probe points

gui_data.bezier_control_points{curr_probe} = ...
    [gui_data.bezier_control_points{curr_probe}; curr_point.Position, gui_data.curr_slice];
[~, sortIdx] = sort(gui_data.bezier_control_points{curr_probe}(:, 3)); % sort points by slice # 1+
gui_data.bezier_control_points{curr_probe} = ...
    gui_data.bezier_control_points{curr_probe}(sortIdx, :);

set(gui_data.histology_ax_title, 'String', ...
    ['Arrows to move, Number to draw probe [', num2str(1), ':', num2str(gui_data.n_probes), '], Esc to save/quit']);

% Plot Bezier curve if there are enough control points (at least 4 for a cubic Bezier curve)
if size(gui_data.bezier_control_points{curr_probe}, 1) >= 4
    t = linspace(0, 1, 100);
    B = bezier_curve(t, gui_data.bezier_control_points{curr_probe});

    % delete any previous Bezier curves
    delete(gui_data.bezier_curves(curr_probe));
    % Compute the new Bezier curve in 3D
    t = linspace(0, 1, 1000);
    B = bezier_curve(t, gui_data.bezier_control_points{gui_data.curr_probe});

    % Filter based on z-slice
    z_slice_tolerance = 0.1; % define a tolerance for how close the z-value of the Bezier curve has to be to z_slice to be displayed
    indices = find(abs(B(:, 3)-gui_data.curr_slice) < z_slice_tolerance); % find indices of points on the Bezier curve close to z_slice
    B_slice = B(indices, :);


    gui_data.bezier_curves(curr_probe) = plot(B_slice(:, 1), B_slice(:, 2), 'Color', gui_data.probe_color(curr_probe, :), 'Parent', gui_data.histology_ax, 'LineWidth', 2);
end

% delete any previous control points and replot
%11delete(gui_data.probe_points{curr_probe});
% gui_data.probe_points{curr_probe} = plot(gui_data.bezier_control_points{curr_probe}(:,1),...
%     gui_data.bezier_control_points{curr_probe}(:,2), 'o', 'MarkerFaceColor', gui_data.probe_color(curr_probe, :),...
%     'MarkerEdgeColor', gui_data.probe_color(curr_probe, :), 'Parent', gui_data.histology_ax);
curr_point.delete;
thesePoints = gui_data.bezier_control_points{curr_probe}(:, 3) == gui_data.curr_slice;
set(gui_data.probe_points{curr_probe}, 'XData', gui_data.bezier_control_points{curr_probe}(thesePoints, 1), 'YData', gui_data.bezier_control_points{curr_probe}(thesePoints, 2), ...
    'MarkerFaceColor', gui_data.probe_color(curr_probe, :), 'MarkerEdgeColor', gui_data.probe_color(curr_probe, :));

% store points
gui_data.probe_points_histology{gui_data.curr_slice, gui_data.curr_probe} = gui_data.bezier_control_points{gui_data.curr_probe};
% Store the gui_data again
guidata(gui_fig, gui_data);

populate_points_table(gui_fig);
end

function update_slice(gui_fig)
% Draw histology and CCF slice

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h, 'CData', (gui_data.slice_im{gui_data.curr_slice})*gui_data.contrast_alpha+gui_data.brightness_beta)

% Clear any current lines, draw probe lines
for iProbe = 1:gui_data.n_probes
    gui_data.bezier_curves(iProbe).delete;
end
%gui_data.probe_fit_lines.delete;
%gui_data.inflection_points_scatter.delete;

for curr_probe = 1:size(gui_data.probe_points_histology, 2)
    if sum(gui_data.bezier_control_points{curr_probe}(:,3) == gui_data.curr_slice) > 0
        thesePoints = gui_data.bezier_control_points{curr_probe}(:, 3) == gui_data.curr_slice;
        set(gui_data.probe_points{curr_probe}, 'XData', gui_data.bezier_control_points{curr_probe}(thesePoints, 1), 'YData', gui_data.bezier_control_points{curr_probe}(thesePoints, 2), ...
            'MarkerFaceColor', gui_data.probe_color(curr_probe, :), 'MarkerEdgeColor', gui_data.probe_color(curr_probe, :));
        t = linspace(0, 1, 1000);
        B = bezier_curve(t, gui_data.bezier_control_points{curr_probe});

        % Filter based on z-slice
        z_slice_tolerance = 0.1; % define a tolerance for how close the z-value of the Bezier curve has to be to z_slice to be displayed
        indices = find(abs(B(:, 3)-gui_data.curr_slice) < z_slice_tolerance); % find indices of points on the Bezier curve close to z_slice
        B_slice = B(indices, :);


        gui_data.bezier_curves(curr_probe) = plot(B_slice(:, 1), B_slice(:, 2), 'Color', gui_data.probe_color(curr_probe, :), 'Parent', gui_data.histology_ax, 'LineWidth', 2);

    else
        try
            set(gui_data.probe_points{curr_probe}, 'XData', NaN, 'YData', NaN);
            set(gui_data.bezier_curves(curr_probe), 'XData', NaN, 'YData', NaN);
        catch
        end


    end
end


set(gui_data.histology_ax_title, 'String', ...
    ['Arrows to move, Number to draw probe [', num2str(1), ':', num2str(gui_data.n_probes), '], Esc to save/quit']);

% Upload gui data
guidata(gui_fig, gui_data);
populate_points_table(gui_fig);

end

function plot_probe(gui_data, probe_ccf)

% Plot probe trajectories
figure('Name', 'Probe trajectories');
axes_atlas = axes;
[~, brain_outline] = plotBrainGrid([], axes_atlas);
set(axes_atlas, 'ZDir', 'reverse');
hold(axes_atlas, 'on');
axis vis3d equal off manual
view([-30, 25]);
caxis([0, 300]);
[ap_max, dv_max, ml_max] = size(gui_data.tv);
% xlim([-10,ml_max+10])
% ylim([-10,dv_max+10])
% zlim([-10,ap_max+10])
h = rotate3d(gca);
h.Enable = 'on';

for curr_probe = 1:length(probe_ccf)
    % Plot points and line of best fit
    slice_points = find(~cellfun(@isempty, gui_data.probe_points_histology(:, curr_probe)'));
    if ~isempty(slice_points)
        linear_fit = 0; %hard coded for now
        spline_fit = 0;
        if linear_fit
            thesePoints = probe_ccf(curr_probe).points * 2.5; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid
            r0 = mean(thesePoints, 1);
            xyz = bsxfun(@minus, thesePoints, r0);
            [~, ~, V] = svd(xyz, 0);
            %V= permute(V, [3, 2, 1]);
            histology_probe_direction = V(:, 1);
            % (make sure the direction goes down in DV - flip if it's going up)
            if histology_probe_direction(2) < 0
                histology_probe_direction = -histology_probe_direction;
            end

            line_eval = [-1000, 1000];
            probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0);
            plot3(thesePoints(:, 1), ...
                thesePoints(:, 2), ...
                thesePoints(:, 3), ...
                '.', 'color', gui_data.probe_color(curr_probe, :), 'MarkerSize', 20);
            line(probe_fit_line(:, 1), probe_fit_line(:, 2), probe_fit_line(:, 3), ...
                'color', gui_data.probe_color(curr_probe, :), 'linewidth', 2)
        elseif spline_fit
            thesePoints = probe_ccf(curr_probe).points * 2.5; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid

            plot3(thesePoints(:, 1), ...
                thesePoints(:, 2), ...
                thesePoints(:, 3), ...
                '.', 'color', gui_data.probe_color(curr_probe, :), 'MarkerSize', 20);

            x = thesePoints(:, 1);
            y = thesePoints(:, 2);
            z = thesePoints(:, 3);

            % (spline)
            t = (1:length(x))'; % Creating a parameter t based on the number of data points

            % Creating fittype objects for smooth spline with smoothing parameter
            ft = fittype('smoothingspline');

            % Define a smoothing parameter
            smoothing_param = gui_data.spline_smoothness_factor(curr_probe); % Adjust this value between 0 and 1 to control the smoothness

            % Fitting the curves with the smoothing parameter
            fitresult_x = fit(t, x, ft, 'SmoothingParam', smoothing_param);
            fitresult_y = fit(t, y, ft, 'SmoothingParam', smoothing_param);
            fitresult_z = fit(t, z, ft, 'SmoothingParam', smoothing_param);

            delta = 5;
            % Step 3: Extending the Curve and Visualizing
            % Extend the range of t for extrapolation
            t_extended = [linspace(min(t)-delta, max(t)+delta, 1000)]'; % Adjust delta to control the extension amount

            % Evaluate the extended spline at the new t values
            x_new = fitresult_x(t_extended);
            y_new = fitresult_y(t_extended);
            z_new = fitresult_z(t_extended);


            plot3(x_new, y_new, z_new, ...
                'color', gui_data.probe_color(curr_probe, :), 'linewidth', 2)
        else
             thesePoints = probe_ccf(curr_probe).points * 2.5; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid

            plot3(thesePoints(:, 1), ...
                thesePoints(:, 2), ...
                thesePoints(:, 3), ...
                '.', 'color', gui_data.probe_color(curr_probe, :), 'MarkerSize', 20);

            x = thesePoints(:, 1);
            y = thesePoints(:, 2);
            z = thesePoints(:, 3);

            % (bezier)
            t = linspace(0, 1, 1000); % Creating a parameter t based on the number of data points
             Bezier_fit = bezier_curve(t,[x, y, z]);
            
              plot3(Bezier_fit(:,1), Bezier_fit(:,2), Bezier_fit(:,3), ...
                'color', gui_data.probe_color(curr_probe, :), 'linewidth', 2)
            

        end

    end
end

% Plot probe areas
figure('Name', 'Trajectory areas');
% (load the colormap - located in the repository, find by associated fcn)
allenCCF_path = fileparts(which('allenCCFbregma'));
cmap_filename = [allenCCF_path, filesep, 'allen_ccf_colormap_2017.mat'];
load(cmap_filename);
for curr_probe = 1:length(probe_ccf)
    slice_points = find(~cellfun(@isempty, gui_data.probe_points_histology(:, curr_probe)'));
    if ~isempty(slice_points)

        curr_axes = subplot(1, gui_data.n_probes, curr_probe);
        % flip probe direction if necessary - check dv (dim 2) goes down
        if probe_ccf(curr_probe).trajectory_coords(1, 2) > probe_ccf(curr_probe).trajectory_coords(end, 2)
            probe_ccf(curr_probe).trajectory_coords = flipud(probe_ccf(curr_probe).trajectory_coords);
            probe_ccf(curr_probe).trajectory_areas = flipud(probe_ccf(curr_probe).trajectory_areas);
            probe_ccf(curr_probe).points = flipud(probe_ccf(curr_probe).points);
        end
        trajectory_area_boundaries = ...
            [1; find(diff(probe_ccf(curr_probe).trajectory_areas) ~= 0); length(probe_ccf(curr_probe).trajectory_areas)];
        trajectory_area_centers = trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries) / 2;

        for iArea = 1:size(trajectory_area_centers, 1)
            trajectory_area_labels(iArea) = gui_data.st.acronym(gui_data.st.id == ...
                probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers(iArea))));
        end
        image(probe_ccf(curr_probe).trajectory_areas);
        colormap(curr_axes, cmap);
        caxis([1, size(cmap, 1)])
        set(curr_axes, 'YTick', trajectory_area_centers, 'YTickLabels', trajectory_area_labels);
        set(curr_axes, 'XTick', []);

        if curr_probe == 1
            title(['Probe ', num2str(curr_probe)]);
        else
            title([num2str(curr_probe)]);
        end
    end

end

end

function updateProbe(gui_fig)
gui_data = guidata(gui_fig);

% Get the selected probe index
selected_probe_idx = get(gui_data.probe_dropdown, 'Value');

% Set the current probe in gui_data
gui_data.curr_probe = selected_probe_idx;

% Refresh the points list
gui_data.selected_row = 1;
gui_data.points_table.Value = gui_data.selected_row;

% Save gui_data
guidata(gui_fig, gui_data);
populate_points_table(gui_fig);

end

function updateProbe_slider(gui_fig, curr_probe)
gui_data = guidata(gui_fig);

% Get the selected probe index
set(gui_data.probe_dropdown, 'Value', curr_probe);


% Refresh the points list
gui_data.selected_row = 1;
gui_data.points_table.Value = gui_data.selected_row;

% Save gui_data
guidata(gui_fig, gui_data);
populate_points_table(gui_fig);


end

function startPtButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Set start button for current probe
curr_point = drawpoint;
gui_data.start_points(gui_data.curr_probe,:) = [curr_point.Position, gui_data.curr_slice];
%gui_data.start_points_scatter{gui_data.curr_probe}.delete;
%gui_data.start_points_scatter{gui_data.curr_probe} = plot(curr_point.Position(1), curr_point.Position(2), 'Color', gui_data.probe_color(gui_data.curr_probe, :))
set(gui_data.start_points_scatter{gui_data.curr_probe}, 'XData',curr_point.Position(1), 'YData',curr_point.Position(2),'Color', gui_data.probe_color(gui_data.curr_probe, :), ...
    'LineWidth', 2)
curr_point.delete;
% Upload gui data
guidata(gui_fig, gui_data);
end

function stopPtButtonPushed(gui_fig)
% Set start button for current probe
curr_point = drawpoint;
gui_data.stop_points(gui_data.curr_probe,:) = [curr_point.Position, gui_data.curr_slice];
%gui_data.start_points_scatter{gui_data.curr_probe}.delete;
%gui_data.start_points_scatter{gui_data.curr_probe} = plot(curr_point.Position(1), curr_point.Position(2), 'Color', gui_data.probe_color(gui_data.curr_probe, :))
set(gui_data.stop_points_scatter{gui_data.curr_probe}, 'XData',curr_point.Position(1), 'YData',curr_point.Position(2),'Color', gui_data.probe_color(gui_data.curr_probe, :), ...
    'LineWidth', 2)
curr_point.delete;
% Upload gui data
guidata(gui_fig, gui_data);
end

%% bezier functions
function selectPoint(gui_fig, table, event)
gui_data = guidata(gui_fig);

removeHighlighted(gui_fig);
%set(gui_data.hHighlighted, 'PickableParts', 'visible');
gui_data.selected_row = event.Source.Value; % Assuming only one row is selectable at a time
selected_point = gui_data.bezier_control_points{gui_data.curr_probe}(gui_data.selected_row, :);
%selected_point = gui_data.probe_points_histology{gui_data.curr_slice, gui_data.curr_probe}{1}(selected_row,:);

% Make the point draggable and set a listener to update the Bezier fit when dragged
set(gui_data.hHighlighted, 'XData', selected_point(:, 1), 'YData', selected_point(:, 2));
gui_data.pointSelected = true; % Indicate that a point is currently selected
%addlistener(selected_point, 'ROIMoved', @(src,evt) bezier_fit_probe(gui_fig));

% Optionally, change the color or appearance of the selected point
%selected_point.Color = 'g';

guidata(gui_fig, gui_data);
end

function removeHighlighted(gui_fig)
gui_data = guidata(gui_fig);
if ishandle(gui_data.hHighlighted)
    set(gui_data.hHighlighted, 'XData', NaN, 'YData', NaN);
end
guidata(gui_fig, gui_data);
end

function populate_points_table(gui_fig)
gui_data = guidata(gui_fig);

%probe_points = gui_data.probe_points_histology{gui_data.curr_slice, gui_data.curr_probe}{1};
%non_empty_slices = find(~cellfun(@isempty, gui_data.probe_points_histology(:, gui_data.curr_probe)));

pointStrings = cell(size(gui_data.bezier_control_points{gui_data.curr_probe}, 1), 1);

if ~isempty(gui_data.bezier_control_points{gui_data.curr_probe})
    for i = 1:size(gui_data.bezier_control_points{gui_data.curr_probe}, 1)
        prefix = '';
        if gui_data.bezier_control_points{gui_data.curr_probe}(i, 3) == gui_data.curr_slice %condition to check if the point is in the current slice%
            prefix = '[CURR] ';
        end
        pointStrings{i} = [prefix, 'Point ', num2str(i), ': (', ...
            num2str(gui_data.bezier_control_points{gui_data.curr_probe}(i, 1)), ', ', num2str(gui_data.bezier_control_points{gui_data.curr_probe}(i, 2)), ')'];
    end
end

set(gui_data.points_table, 'String', pointStrings);
guidata(gui_fig, gui_data);
end


function deleteSelectedPoint(~, ~, gui_fig)
gui_data = guidata(gui_fig);
selectedIdx = gui_data.selected_row;
if ~isempty(selectedIdx)
    gui_data.bezier_control_points{gui_data.curr_probe}(selectedIdx, :) = [];
    gui_data.probe_points_histology{gui_data.curr_slice, gui_data.curr_probe}(selectedIdx, :) = [];
end
gui_data.pointSelected = true;
guidata(gui_fig, gui_data);
removeHighlighted(gui_fig);
if gui_data.selected_row > 1
    gui_data.selected_row = gui_data.selected_row - 1;
else
    gui_data.selected_row = gui_data.selected_row + 1;
end
gui_data.points_table.Value = gui_data.selected_row;
guidata(gui_fig, gui_data);
populate_points_table(gui_fig)
update_slice(gui_fig);

end

function startDragFcn(~, ~, gui_fig)
gui_data = guidata(gui_fig);
selectedIdx = gui_data.selected_row;
if ~isempty(selectedIdx)
    set(gui_fig, 'WindowButtonMotionFcn', @(src, event)draggingFcn(src, event, gui_fig));
    set(gui_fig, 'WindowButtonUpFcn', @(src, event)stopDragFcn(src, event, gui_fig));
end
end

function draggingFcn(~, ~, gui_fig)
gui_data = guidata(gui_fig);

% Assuming your draggable points are on the current axis (gca).
% If they are on a specific axis, replace gca with the handle to that axis.
currentAxis = gca;

% Get the current point with respect to the intended axis
pt = get(currentAxis, 'CurrentPoint');
draggedPos = pt(1, 1:2);

% Get the axis limits
xLim = get(currentAxis, 'XLim');
yLim = get(currentAxis, 'YLim');

% Check and adjust the x position if it's outside the axis limits
if draggedPos(1) < xLim(1)
    draggedPos(1) = xLim(1);
elseif draggedPos(1) > xLim(2)
    draggedPos(1) = xLim(2);
end

% Check and adjust the y position if it's outside the axis limits
if draggedPos(2) < yLim(1)
    draggedPos(2) = yLim(1);
elseif draggedPos(2) > yLim(2)
    draggedPos(2) = yLim(2);
end

% Update the position of the point
selectedIdx = gui_data.selected_row;
gui_data.bezier_control_points{gui_data.curr_probe}(selectedIdx, 1:2) = draggedPos;
gui_data.probe_points_histology{gui_data.curr_slice, gui_data.curr_probe}(selectedIdx, 1:2) = draggedPos;

set(gui_data.hHighlighted, 'XData', draggedPos(1), 'YData', draggedPos(2));
gui_data.bezier_curves(gui_data.curr_probe).delete;

t = linspace(0, 1, 1000);
B = bezier_curve(t, gui_data.bezier_control_points{gui_data.curr_probe});

% Filter based on z-slice
z_slice_tolerance = 0.1; % define a tolerance for how close the z-value of the Bezier curve has to be to z_slice to be displayed
indices = find(abs(B(:, 3)-gui_data.curr_slice) < z_slice_tolerance); % find indices of points on the Bezier curve close to z_slice
B_slice = B(indices, :);


gui_data.bezier_curves(gui_data.curr_probe) = plot(B_slice(:, 1), B_slice(:, 2), 'Color', gui_data.probe_color(gui_data.curr_probe, :), 'Parent', gui_data.histology_ax, 'LineWidth', 2);


guidata(gui_fig, gui_data);
populate_points_table(gui_fig);
removeHighlighted(gui_fig);

end

function stopDragFcn(~, ~, gui_fig)
gui_data = guidata(gui_fig);
set(gui_fig, 'WindowButtonMotionFcn', '');
set(gui_fig, 'WindowButtonUpFcn', '');
gui_data.pointSelected = false; % Reset the flag since we're no longer selecting a point
guidata(gui_fig, gui_data);
update_slice(gui_fig);
end


% Bezier curve function
function B = bezier_curve_prev(t, control_points)
n = size(control_points, 1) - 1; % degree of the polynomial
B = zeros(2, length(t));

for i = 0:n
    B = B + nchoosek(n, i) * (1 - t).^(n - i) .* t.^i .* control_points(i+1, :)';
end
B = B';
end

function B = bezier_curve(t, control_points)
n = size(control_points, 1) - 1; % degree of the polynomial
B = zeros(3, length(t)); % Change to 3 for 3D

[~, control_points_idx] = sort(control_points(:,3));
control_points = control_points(control_points_idx, :);
for i = 0:n
    B = B + nchoosek(n, i) * (1 - t).^(n - i) .* t.^i .* control_points(i+1, :)';
end
B = B';
end

%% control appearance

function autoContrastButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Set contrast and brightness
% Auto based on this function :
% g(i,j) = α * f(i,j) + β - we want to find α (contrast) and β
% (brightness)
% α = 255 / (maximum_gray - minimum_gray)
% solve β by plugging it into the formula where g(i, j)=0 and
% f(i, j)=minimum_gray

gui_data.contrast_alpha = 255 / double(max(gui_data.slice_im{gui_data.curr_slice}, [], 'all')- ...
    min(gui_data.slice_im{gui_data.curr_slice}, [], 'all'));
gui_data.brightness_beta = -(gui_data.contrast_alpha * double(min(gui_data.slice_im{gui_data.curr_slice}, [], 'all')));

% update bright and contast slider values
gui_data.brightness_slider.Value = gui_data.brightness_beta;
gui_data.contrast_slider.Value = gui_data.contrast_alpha;

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_slice(gui_fig)

end

function brightnessButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

brightnessValue = gui_data.brightness_slider.Value;

gui_data.brightness_beta = brightnessValue;

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_slice(gui_fig)

end

function pointSizeButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

pointSize = gui_data.point_size_slider.Value;

gui_data.point_size = pointSize;

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_slice(gui_fig)

end

function contrastButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

contrastValue = gui_data.contrast_slider.Value;

gui_data.contrast_alpha = contrastValue;

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_slice(gui_fig)

end

%% main buttons
function saveButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Save current data
% Initialize structure to save
probe_ccf = struct( ...
    'points', cell(gui_data.n_probes, 1), ...
    'trajectory_coords', cell(gui_data.n_probes, 1), ... .
    'trajectory_areas', cell(gui_data.n_probes, 1));

%ccf_slice_fn = ['/home/netshare/zaru/JF096/Histology/downsampled_stacks/025_micron/brainReg/manual/histology_ccf.mat'];
%load(ccf_slice_fn);
%gui_data.histology_ccf = histology_ccf;

% Convert probe points to CCF points by alignment and save
for curr_probe = 1:gui_data.n_probes
    slice_points = find(~cellfun(@isempty, gui_data.probe_points_histology(:, curr_probe)'));
    if ~isempty(slice_points)
        for curr_slice = slice_points


            probe_points_atlas_x = gui_data.probe_points_histology{curr_slice, curr_probe}(:, 1);
            probe_points_atlas_y = gui_data.probe_points_histology{curr_slice, curr_probe}(:, 2);

            probe_points_atlas_x = round(probe_points_atlas_x);
            probe_points_atlas_y = round(probe_points_atlas_y);

            % Get CCF coordinates corresponding to atlas slice points
            % (CCF coordinates are in [AP,DV,ML])
            use_points = find(~isnan(probe_points_atlas_x) & ~isnan(probe_points_atlas_y));


            for curr_point = 1:length(use_points)
                ccf_ap = gui_data.histology_ccf(curr_slice). ...
                    plane_ap(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                ccf_ml = gui_data.histology_ccf(curr_slice). ...
                    plane_dv(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                ccf_dv = gui_data.histology_ccf(curr_slice). ...
                    plane_ml(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                probe_ccf(curr_probe).points = ...
                    vertcat(probe_ccf(curr_probe).points, [ccf_ap, ccf_dv, ccf_ml]);
            end
        end

        % Sort probe points by DV (probe always top->bottom)
        [~, dv_sort_idx] = sort(probe_ccf(curr_probe).points(:, 2));
        probe_ccf(curr_probe).points = probe_ccf(curr_probe).points(dv_sort_idx, :);
    end

end
gui_data.probe_ccf = probe_ccf;
% Get areas along probe trajectory

for curr_probe = 1:gui_data.n_probes
    % Get slice and point info
    slice_points = find(~cellfun(@isempty, gui_data.probe_points_histology(:, curr_probe)')); % slices that contain points

    if ~isempty(slice_points)
        for curr_slice = slice_points


            probe_points_atlas_x = gui_data.probe_points_histology{curr_slice, curr_probe}(:, 1);
            probe_points_atlas_y = gui_data.probe_points_histology{curr_slice, curr_probe}(:, 2);

            probe_points_atlas_x = round(probe_points_atlas_x);
            probe_points_atlas_y = round(probe_points_atlas_y);

            % Get CCF coordinates corresponding to atlas slice points
            % (CCF coordinates are in [AP,DV,ML])
            use_points = find(~isnan(probe_points_atlas_x) & ~isnan(probe_points_atlas_y));
            for curr_point = 1:length(use_points)
                ccf_ap = gui_data.histology_ccf(curr_slice). ...
                    plane_ap(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                ccf_ml = gui_data.histology_ccf(curr_slice). ...
                    plane_dv(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                ccf_dv = gui_data.histology_ccf(curr_slice). ...
                    plane_ml(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                gui_data.probe_ccf(curr_probe).points = ...
                    vertcat(gui_data.probe_ccf(curr_probe).points, [ccf_ap, ccf_dv, ccf_ml]);
            end
        end

        % Sort probe points by DV (probe always top->bottom)
        [~, dv_sort_idx] = sort(gui_data.probe_ccf(curr_probe).points(:, 2));
        gui_data.probe_ccf(curr_probe).points = gui_data.probe_ccf(curr_probe).points(dv_sort_idx, :);

        linear_fit = 0; % QQ hard-coded for now
        spline_fit = 0;
        if linear_fit

            % Get best fit line through points as probe trajectory
            r0 = mean(gui_data.probe_ccf(curr_probe).points, 1);
            xyz = bsxfun(@minus, gui_data.probe_ccf(curr_probe).points, r0);
            [~, ~, V] = svd(xyz, 0);
            histology_probe_direction = V(:, 1);

            % (make sure the direction goes down in DV - flip if it's going up)
            if histology_probe_direction(3) < 0
                histology_probe_direction = -histology_probe_direction;
            end

            line_eval = [-1000, 1000];
            probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0)';

            % Get the positions of the probe trajectory
            trajectory_n_coords = max(abs(diff(probe_fit_line, [], 2)));

            % get AP, DV, ML
            [trajectory_ap_ccf, trajectory_dv_ccf, trajectory_ml_ccf] = deal( ...
                round(linspace(probe_fit_line(1, 1), probe_fit_line(1, 2), trajectory_n_coords)), ...
                round(linspace(probe_fit_line(3, 1), probe_fit_line(3, 2), trajectory_n_coords)), ...
                round(linspace(probe_fit_line(2, 1), probe_fit_line(2, 2), trajectory_n_coords)));


        elseif spline_fit % spline fit
            ap_points = gui_data.probe_ccf(curr_probe).points(:, 1);
            dv_points = gui_data.probe_ccf(curr_probe).points(:, 2);
            ml_points = gui_data.probe_ccf(curr_probe).points(:, 3);

            t = (1:length(ap_points))'; % Creating a parameter t based on the number of data points

            % Creating fittype objects for smooth spline with smoothing parameter
            ft = fittype('smoothingspline');

            % Define a smoothing parameter
            smoothing_param = gui_data.spline_smoothness_factor(curr_probe); % Adjust this value between 0 and 1 to control the smoothness

            % Fitting the curves with the smoothing parameter
            fitresult_ap = fit(t, ap_points, ft, 'SmoothingParam', smoothing_param);
            fitresult_dv = fit(t, dv_points, ft, 'SmoothingParam', smoothing_param);
            fitresult_ml = fit(t, ml_points, ft, 'SmoothingParam', smoothing_param);

            delta = 5;
            % (Extending the Curve and Visualizing)

            % Extend the range of t for extrapolation
            t_extended = [linspace(min(t)-delta, max(t)+delta, 1000)]'; % Adjust delta to control the extension amount

            % Evaluate the extended spline at the new t values
            trajectory_ap_ccf = fitresult_ap(t_extended)';
            trajectory_ml_ccf = fitresult_ml(t_extended)';
            trajectory_dv_ccf = fitresult_dv(t_extended)';

        else % manual bezier curve drawing 
            ap_points = gui_data.probe_ccf(curr_probe).points(:, 1);
            dv_points = gui_data.probe_ccf(curr_probe).points(:, 2);
            ml_points = gui_data.probe_ccf(curr_probe).points(:, 3);

            t = linspace(0, 1, 1000);
            Bezier_fit = bezier_curve(t,[ap_points, dv_points, ml_points]);

            trajectory_ap_ccf = Bezier_fit(:,1)';
            trajectory_ml_ccf = Bezier_fit(:,2)';
            trajectory_dv_ccf = Bezier_fit(:,3)';


            

        end

        % remove any out-of-bounds points
        trajectory_coords_outofbounds = ...
            any([trajectory_ap_ccf; trajectory_dv_ccf; trajectory_ml_ccf] < 1, 1) | ...
            any([trajectory_ap_ccf; trajectory_dv_ccf; trajectory_ml_ccf] > size(gui_data.av)', 1);

        trajectory_coords = ...
            [trajectory_ap_ccf(~trajectory_coords_outofbounds)', ...
            trajectory_dv_ccf(~trajectory_coords_outofbounds)', ...
            trajectory_ml_ccf(~trajectory_coords_outofbounds)'];

        trajectory_coords_idx = sub2ind(size(gui_data.av), ...
            round(trajectory_coords(:, 1)), round(trajectory_coords(:, 2)), round(trajectory_coords(:, 3)));

        trajectory_areas_uncut = gui_data.av(trajectory_coords_idx)';

        % Get rid of NaN's and start/end 1's (non-parsed)
        trajectory_areas_parsed = find(trajectory_areas_uncut > 1); %ones(size(trajectory_areas_uncut,2),1);%find(trajectory_areas_uncut > 1);
        use_trajectory_areas = trajectory_areas_parsed(1): ...
            trajectory_areas_parsed(end);
        trajectory_areas = reshape(trajectory_areas_uncut(use_trajectory_areas), [], 1);

        % store fit points in gui_data
        gui_data.probe_ccf(curr_probe).trajectory_coords = double(trajectory_coords(use_trajectory_areas, :));
        gui_data.probe_ccf(curr_probe).trajectory_areas = double(trajectory_areas);

    end


end

% Save fitted probe CCF points
save_fn = [gui_data.slice_im_path, filesep, 'probe_ccf.mat'];
save(save_fn, 'probe_ccf');

% save points
save_fn_pt = [gui_data.slice_im_path, filesep, 'probe_points.mat'];
probe_points = gui_data.probe_points_histology;
save(save_fn_pt, 'probe_points');

% save fit info
% save_fn_fits = [gui_data.slice_im_path, filesep, 'probe_fit_sline_smoothness.mat'];
% probe_fit_sline_smoothness = gui_data.spline_smoothness_factor;
% save(save_fn_fits, 'probe_fit_sline_smoothness');

% save start + stop points 
save_fn_pt = [gui_data.slice_im_path, filesep, 'probe_points_start.mat'];
probe_points_start = gui_data.start_points;
save(save_fn_pt, 'probe_points_start');

save_fn_pt = [gui_data.slice_im_path, filesep, 'probe_points_stop.mat'];
probe_points_stop = gui_data.stop_points;
save(save_fn_pt, 'probe_points_stop');

% Display
disp(['Saved probe locations in ', save_fn])

% Plot probe trajectories
plot_probe(gui_data, gui_data.probe_ccf);

% Upload gui data
guidata(gui_fig, gui_data);

end

function loadButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% load

save_fn = [gui_data.slice_im_path, filesep, 'probe_ccf.mat'];
load(save_fn);

% load points
try
    save_fn_pt = [gui_data.slice_im_path, filesep, 'probe_points.mat'];
    load(save_fn_pt);
catch
    save_fn_pt = [gui_data.slice_im_path, filesep, 'probe_points_histology.mat'];
    load(save_fn_pt);
    probe_points = pp_histo;
end

% load start/stop points if present 

try
    save_fn_pt = [gui_data.slice_im_path, filesep, 'probe_points_start.mat'];
    load(save_fn_pt);
    save_fn_pt = [gui_data.slice_im_path, filesep, 'probe_points_stop.mat'];
    load(save_fn_pt);
    gui_data.stop_points = probe_points_stop;
    gui_data.start_points = probe_points_start;
catch
    gui_data.stop_points = nan(gui_data.n_probes, 3);
    gui_data.start_points = nan(gui_data.n_probes, 3);
end



% slices AP
for iSlice = 1:size(gui_data.histology_ccf, 1)
    ap_all(iSlice) = gui_data.histology_ccf(iSlice).plane_ap(1);
end
previous_string = gui_data.histology_ax_title.String;
set(gui_data.histology_ax_title, 'String', 'loading...')
previous_n_probes = gui_data.n_probes;
% Add probes if needed
if size(probe_ccf, 1) > gui_data.n_probes


    % Add probes
    gui_data.n_probes = size(probe_ccf, 1);
    gui_data.probe_color = lines(size(probe_ccf, 1));
    for iNewProbe = previous_n_probes + 1:size(probe_ccf, 1)
        gui_data.probe_points_histology(:, iNewProbe) = cell(length(gui_data.slice_im), 1);
        gui_data.probe_points{iNewProbe} = plot(NaN, NaN);
        gui_data.probe_fit_lines(iNewProbe) = gobjects(1, 1);
        gui_data.prev_spline_fit(iNewProbe) = 0.1;
        gui_data.bezier_control_points(iNewProbe) = cell(1, 1);
        gui_data.bezier_curves(iNewProbe) = gobjects(1, 1);
        gui_data.probe_dropdown.String{iNewProbe} = ['Probe ' num2str(iNewProbe)];
        gui_data.start_points_scatter{iNewProbe} = plot(NaN, NaN, 'rx');
        gui_data.stop_points_scatter{iNewProbe} = plot(NaN, NaN, 'rx');
    end

    % Update all buttons

    

end


gui_data.probe_points_histology = probe_points;
for iProbe = 1:gui_data.n_probes
    non_empty_slices = find(~cellfun(@isempty, gui_data.probe_points_histology(:, iProbe)));
    for iSlice = 1:size(non_empty_slices,1)
        gui_data.bezier_control_points{iProbe} = [gui_data.bezier_control_points{iProbe}; ...
            probe_points{non_empty_slices(iSlice), iProbe},...
            repmat(non_empty_slices(iSlice), length(probe_points{non_empty_slices(iSlice), iProbe}), 1)];
    end
    [~, idx] = sort(gui_data.bezier_control_points{iProbe}(:,3));
    gui_data.bezier_control_points{iProbe} =gui_data.bezier_control_points{iProbe}(idx, :);
    
end
set(gui_data.histology_ax_title, 'String', 'successfully loaded')
disp('loaded data')

% Upload gui data
set(gui_data.histology_ax_title, 'String', previous_string)
guidata(gui_fig, gui_data);
populate_points_table(gui_fig)

end

function addProbe(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Add a probe
gui_data.n_probes = gui_data.n_probes + 1;
gui_data.probe_color = lines(gui_data.n_probes);
gui_data.probe_points_histology(:, gui_data.n_probes) = cell(length(gui_data.slice_im), 1);
gui_data.probe_points{gui_data.n_probes} = gobjects(10, 1);
gui_data.probe_inflection_pts(gui_data.n_probes) = gobjects(1, 1);
gui_data.probe_fit_lines(gui_data.n_probes) = gobjects(1, 1);
gui_data.prev_spline_fit(gui_data.n_probes) = 0.1;
gui_data.bezier_control_points = [gui_data.bezier_control_points; cell(1, 1)];
gui_data.bezier_curves(gui_data.n_probes) = plot(NaN, NaN);
gui_data.probe_dropdown.String{gui_data.n_probes} = ['Probe ', num2str(gui_data.n_probes)];
gui_data.start_points_scatter{gui_data.n_probes} = plot(NaN, NaN, 'rx');
gui_data.stop_points_scatter{gui_data.n_probes} = plot(NaN, NaN, 'rx');
gui_data.stop_points(gui_data.n_probes,3) = nan(1, 3);
gui_data.start_points(gui_data.n_probes,3) = nan(1, 3);
% Upload gui data
guidata(gui_fig, gui_data);

end

function toggleVisiblityProbeButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% hide/show
if contains(gui_data.del_probe_btns.String, 'ide')
    gui_data.del_probe_btns.String = 'show'; %<HTML><center><FONT color="white"><b>
    gui_data.visibility = 0;
    delete(gui_data.probe_lines(gui_data.curr_probe)); %QQ fix

else
    curr_probe = gui_data.curr_probe;
    gui_data.del_probe_btns.String = 'hide';
    gui_data.visibility = 1;
    for iPoint = 1:size(gui_data.probe_points_histology{gui_data.curr_slice, curr_probe}, 1)
        gui_data.probe_points{curr_probe}(iPoint) = ...
            scatter(gui_data.probe_points_histology{gui_data.curr_slice, curr_probe}(iPoint, 1), ...
            gui_data.probe_points_histology{gui_data.curr_slice, curr_probe}(iPoint, 2), ...
            gui_data.point_size, [gui_data.probe_color(curr_probe, :)], 'filled');
    end
end
% Upload gui data
guidata(gui_fig, gui_data);

end

function toggleAllProbeButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% hide/show
if contains(gui_data.toggle_probe_btn.String, 'Hide')
    gui_data.toggle_probe_btn.String = 'Show all other probes';
    other_probes = logical(ones(gui_data.n_probes, 1));
    other_probes(gui_data.curr_probe) = 0;
    other_probes = find(other_probes);
    gui_data.visibility = 0;
    for iProbe = 1:size(other_probes, 1)
        thisProbe = other_probes(iProbe);
        delete(gui_data.probe_points{thisProbe}); %QQ fix
    end
else
    gui_data.toggle_probe_btn.String = 'Hide all other probes';
    other_probes = logical(ones(gui_data.n_probes, 1));
    other_probes(gui_data.curr_probe) = 0;
    other_probes = find(other_probes);
    gui_data.visibility = 1;
    for iProbe = 1:size(other_probes, 1)
        curr_probe = other_probes(iProbe);
        % update probe line
        for iPoint = 1:size(gui_data.probe_points_histology{gui_data.curr_slice, curr_probe}, 1)
            gui_data.probe_points{curr_probe}(iPoint) = ...
                scatter(gui_data.probe_points_histology{gui_data.curr_slice, curr_probe}(iPoint, 1), ...
                gui_data.probe_points_histology{gui_data.curr_slice, curr_probe}(iPoint, 2), ...
                gui_data.point_size, [gui_data.probe_color(curr_probe, :)], 'filled');
        end
    end
end
% Upload gui data
guidata(gui_fig, gui_data);
end



%% local probe buttons

function viewFitButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Get current probe
curr_probe = find([gui_data.viewFit(:).Value]);

% hide/show
if contains(gui_data.viewFit(curr_probe).String, 'iew')
    % set visibility to view
    gui_data.viewFit(curr_probe).String = 'Hide fit'; % to make white: <HTML><center><FONT color="white"><b>
    gui_data.fit_visibility = 1;


    % Convert probe points to CCF points by alignment
    gui_data.probe_ccf = struct( ...
        'points', cell(gui_data.n_probes, 1), ...
        'trajectory_coords', cell(gui_data.n_probes, 1), ... .
        'trajectory_areas', cell(gui_data.n_probes, 1));

    % Get slice and point info
    slice_points = find(~cellfun(@isempty, gui_data.probe_points_histology(:, curr_probe)')); % slices that contain points

    if ~isempty(slice_points)
        for curr_slice = slice_points


            probe_points_atlas_x = gui_data.probe_points_histology{curr_slice, curr_probe}(:, 1);
            probe_points_atlas_y = gui_data.probe_points_histology{curr_slice, curr_probe}(:, 2);

            probe_points_atlas_x = round(probe_points_atlas_x);
            probe_points_atlas_y = round(probe_points_atlas_y);

            % Get CCF coordinates corresponding to atlas slice points
            % (CCF coordinates are in [AP,DV,ML])
            use_points = find(~isnan(probe_points_atlas_x) & ~isnan(probe_points_atlas_y));
            for curr_point = 1:length(use_points)
                ccf_ap = gui_data.histology_ccf(curr_slice). ...
                    plane_ap(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                ccf_ml = gui_data.histology_ccf(curr_slice). ...
                    plane_dv(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                ccf_dv = gui_data.histology_ccf(curr_slice). ...
                    plane_ml(probe_points_atlas_y(curr_point), ...
                    probe_points_atlas_x(curr_point));
                gui_data.probe_ccf(curr_probe).points = ...
                    vertcat(gui_data.probe_ccf(curr_probe).points, [ccf_ap, ccf_dv, ccf_ml]);
            end
        end

        % Sort probe points by DV (probe always top->bottom)
        [~, dv_sort_idx] = sort(gui_data.probe_ccf(curr_probe).points(:, 2));
        gui_data.probe_ccf(curr_probe).points = gui_data.probe_ccf(curr_probe).points(dv_sort_idx, :);

        linear_fit = 0; % QQ hard-coded for now
        if linear_fit

            % Get best fit line through points as probe trajectory
            r0 = mean(gui_data.probe_ccf(curr_probe).points, 1);
            xyz = bsxfun(@minus, gui_data.probe_ccf(curr_probe).points, r0);
            [~, ~, V] = svd(xyz, 0);
            histology_probe_direction = V(:, 1);

            % (make sure the direction goes down in DV - flip if it's going up)
            if histology_probe_direction(3) < 0
                histology_probe_direction = -histology_probe_direction;
            end

            line_eval = [-1000, 1000];
            probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0)';

            % Get the positions of the probe trajectory
            trajectory_n_coords = max(abs(diff(probe_fit_line, [], 2)));

            % get AP, DV, ML
            [trajectory_ap_ccf, trajectory_dv_ccf, trajectory_ml_ccf] = deal( ...
                round(linspace(probe_fit_line(1, 1), probe_fit_line(1, 2), trajectory_n_coords)), ...
                round(linspace(probe_fit_line(3, 1), probe_fit_line(3, 2), trajectory_n_coords)), ...
                round(linspace(probe_fit_line(2, 1), probe_fit_line(2, 2), trajectory_n_coords)));


        else % spline fit
            ap_points = gui_data.probe_ccf(curr_probe).points(:, 1);
            dv_points = gui_data.probe_ccf(curr_probe).points(:, 2);
            ml_points = gui_data.probe_ccf(curr_probe).points(:, 3);

            t = (1:length(ap_points))'; % Creating a parameter t based on the number of data points

            % Creating fittype objects for smooth spline with smoothing parameter
            ft = fittype('smoothingspline');

            % Define a smoothing parameter
            smoothing_param = gui_data.spline_smoothness_factor(curr_probe); % Adjust this value between 0 and 1 to control the smoothness

            % Fitting the curves with the smoothing parameter
            fitresult_ap = fit(t, ap_points, ft, 'SmoothingParam', smoothing_param);
            fitresult_dv = fit(t, dv_points, ft, 'SmoothingParam', smoothing_param);
            fitresult_ml = fit(t, ml_points, ft, 'SmoothingParam', smoothing_param);

            delta = 5;
            % (Extending the Curve and Visualizing)

            % Extend the range of t for extrapolation
            t_extended = [linspace(min(t)-delta, max(t)+delta, 1000)]'; % Adjust delta to control the extension amount

            % Evaluate the extended spline at the new t values
            trajectory_ap_ccf = fitresult_ap(t_extended)';
            trajectory_dv_ccf = fitresult_ml(t_extended)';
            trajectory_ml_ccf = fitresult_dv(t_extended)';

        end

        % remove any out-of-bounds points
        trajectory_coords_outofbounds = ...
            any([trajectory_ap_ccf; trajectory_dv_ccf; trajectory_ml_ccf] < 1, 1) | ...
            any([trajectory_ap_ccf; trajectory_dv_ccf; trajectory_ml_ccf] > size(gui_data.av)', 1);

        trajectory_coords = ...
            [trajectory_ap_ccf(~trajectory_coords_outofbounds)', ...
            trajectory_dv_ccf(~trajectory_coords_outofbounds)', ...
            trajectory_ml_ccf(~trajectory_coords_outofbounds)'];

        trajectory_coords_idx = round(sub2ind(size(gui_data.av), ...
            trajectory_coords(:, 1), trajectory_coords(:, 2), trajectory_coords(:, 3)));

        trajectory_areas_uncut = gui_data.av(trajectory_coords_idx)';

        % Get rid of NaN's and start/end 1's (non-parsed)
        trajectory_areas_parsed = find(trajectory_areas_uncut > 1); %ones(size(trajectory_areas_uncut,2),1);%find(trajectory_areas_uncut > 1);
        use_trajectory_areas = trajectory_areas_parsed(1): ...
            trajectory_areas_parsed(end);
        trajectory_areas = reshape(trajectory_areas_uncut(use_trajectory_areas), [], 1);

        % store fit points in gui_data
        gui_data.probe_ccf(curr_probe).trajectory_coords = double(trajectory_coords(use_trajectory_areas, :));
        gui_data.probe_ccf(curr_probe).trajectory_areas = double(trajectory_areas);

    end

    % find line on this slice
    slice_coords_fit = gui_data.probe_ccf(curr_probe).trajectory_coords( ...
        round(gui_data.probe_ccf(curr_probe).trajectory_coords(:, 1)) == gui_data.curr_slice, 2:3);
    if ~isempty(round(gui_data.probe_ccf(curr_probe).trajectory_coords(:, 1)) == gui_data.curr_slice)
        try
            gui_data.probe_fit_lines(curr_probe) = line([slice_coords_fit(1, 2), slice_coords_fit(end, 2)], ...
                [slice_coords_fit(1, 1), slice_coords_fit(end, 1)], ...
                'linewidth', 3, 'LineStyle', '-', 'color', [gui_data.probe_color(curr_probe, :), 0.5]);
        catch
        end
    end

else
    gui_data.viewFit(curr_probe).String = 'view fit';
    gui_data.fit_visibility = 0;
    try
        gui_data.probe_fit_lines(curr_probe).Color(4) = 0;
    catch
    end

end

% Upload gui data
guidata(gui_fig, gui_data);
end

function splineSmoothnessButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Get current probe
curr_probe = find(diff([gui_data.prev_spline_fit; gui_data.spline_smoothness(:).Value]));

% Update spline smoothness
spineSmoothValue = gui_data.spline_smoothness(curr_probe).Value;
gui_data.prev_spline_fit(curr_probe) = gui_data.spline_smoothness(curr_probe).Value;
gui_data.spline_smoothness(curr_probe).Value = spineSmoothValue;

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_slice(gui_fig)

end

function resetProbeSliceButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Get current probe
curr_probe = find([gui_data.reset_probe_btns(:).Value]);

% Delete current probe histology points and line
gui_data.probe_points_histology{gui_data.curr_slice, curr_probe} = [];
delete(gui_data.probe_points{curr_probe})
%gui_data.probe_points{curr_probe} = gobjects(10,1);
delete(gui_data.probe_fit_lines(curr_probe))
% Upload gui data
guidata(gui_fig, gui_data);

end

function resetProbeGlobalButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Get current probe
curr_probe = find([gui_data.reset_all_probe_btns(:).Value]);

% Delete current probe histology points and line
for iSlice = 1:length(gui_data.slice_im)
    gui_data.probe_points_histology{iSlice, curr_probe} = [];
end
delete(gui_data.probe_points{curr_probe})
delete(gui_data.probe_fit_lines(curr_probe))
% Upload gui data
guidata(gui_fig, gui_data);

end

function resetGlobalButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);


% Delete current probe histology points and line
for iSlice = 1:length(gui_data.slice_im)
    for iProbe = 1:gui_data.n_probes
        gui_data.probe_points_histology{iSlice, iProbe} = [];
    end
end
for iProbe = 1:gui_data.n_probes
    delete(gui_data.probe_points{iProbe})
end
delete(gui_data.probe_fit_lines(:))

% Upload gui data
guidata(gui_fig, gui_data);

end
