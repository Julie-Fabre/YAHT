function ya_alignEphysAndHistology_draggable(st, outputDir, ...
    spike_times, spike_templates, template_depths, spike_xdepths, template_xdepths, ...
    lfp, lfp_channel_positions, lfp_channel_xpositions, use_probe, probeLength, curr_shank, unitType)
% based on the great AP_align_probe_histology(st,slice_path,spike_times,spike_templates,template_depths,lfp,lfp_channel_positions,use_probe)

% If no probe specified, use probe 1
if ~exist('use_probe', 'var') || isempty(use_probe)
    use_probe = 1;
end

% Load probe CCF
probe_ccf_fn = [outputDir, filesep, 'probe_ccf.mat'];
load(probe_ccf_fn);

if ~isnan(curr_shank)
    theseChannelPositions = [(curr_shank - 1) * 250, (curr_shank - 1) * 250 + 32];
    theseChannels = ismember(lfp_channel_xpositions, theseChannelPositions);
    theseTemplates = ismember(template_xdepths, theseChannelPositions);
    theseSpikes = ismember(spike_xdepths, theseChannelPositions);
    spike_times = spike_times(theseSpikes);
    spike_templates = spike_templates(theseSpikes);
    %rename
    good_templates_idx = unique(spike_templates);
    new_spike_idx = nan(max(spike_templates), 1);
    new_spike_idx(good_templates_idx) = 1:length(good_templates_idx);
    spike_templates = new_spike_idx(spike_templates);
    unitType = unitType(theseTemplates);
    template_depths = template_depths(theseTemplates);
    [~, ~, spike_templates_reidx] = unique(spike_templates);
    norm_template_spike_n = mat2gray(log10(accumarray(spike_templates_reidx, 1)+1));

else
    theseChannels = ones(size(lfp_channel_xpositions, 1), 1);
    theseTemplates = ones(size(template_xdepths, 1), 1);
    theseSpikes = ones(size(spike_xdepths, 1), 1);
    [~, ~, spike_templates_reidx] = unique(spike_templates);
    norm_template_spike_n = mat2gray(log10(accumarray(spike_templates_reidx, 1)+1));
end
% Get normalized log spike n


% Get multiunit correlation
if probeLength < 1200
    n_corr_groups = 20;
else
    n_corr_groups = 40;
end
max_depths = probeLength;
% if min(template_depths) < 1440 % single shank 
    min_depths = 0;
% elseif min(template_depths) < 2175 % 2 shanks
%     min_depths = 1440;
% else % 4shanks
%     min_depths = 2175;
% end
depth_group_edges = linspace(min_depths, max_depths, n_corr_groups+1);


depth_group = discretize(template_depths, depth_group_edges);
depth_group_centers = depth_group_edges(1:end-1) + (diff(depth_group_edges) / 2);
unique_depths = 1:length(depth_group_edges) - 1;

spike_binning = 0.01; % seconds
corr_edges = nanmin(spike_times):spike_binning:nanmax(spike_times);
corr_centers = corr_edges(1:end-1) + diff(corr_edges);

binned_spikes_depth = zeros(length(unique_depths), length(corr_edges)-1);
for curr_depth = 1:length(unique_depths)
    binned_spikes_depth(curr_depth, :) = histcounts(spike_times( ...
        ismember(spike_templates, find(depth_group == unique_depths(curr_depth)))), ...
        corr_edges);
end

mua_corr = corrcoef(binned_spikes_depth');
gui_data = struct;
gui_fig = figure('color', 'w', 'KeyPressFcn', @keypress);

% Plot spike depth vs rate
unit_ax = subplot('Position', [0.1, 0.1, 0.1, 0.8]);
scatter(norm_template_spike_n, template_depths, 15, 'k', 'filled');
set(unit_ax, 'YDir', 'reverse');
ylim([min_depths, max_depths]);
xlabel('N spikes')
title('Template depth & rate')
set(unit_ax, 'FontSize', 12)
ylabel('Depth (\mum)');

% Plot multiunit correlation
multiunit_ax = subplot('Position', [0.2, 0.1, 0.2, 0.8]);
imagesc(depth_group_centers, depth_group_centers, mua_corr);
%caxis([0,max(mua_corr(mua_corr ~= 1))]); colormap(hot);
ylim([min_depths, max_depths]);
set(multiunit_ax, 'YTick', []);
title('MUA correlation');
set(multiunit_ax, 'FontSize', 12)
xlabel(multiunit_ax, 'Multiunit depth');

binSize = n_corr_groups;
try
units_perDepth = histcounts(template_depths, min(template_depths):binSize:max(template_depths));
noiseUnits_perDepth = histcounts(template_depths(unitType == 0), min(template_depths):binSize:max(template_depths))./units_perDepth;
goodUnits_perDepth = histcounts(template_depths(unitType == 1), min(template_depths):binSize:max(template_depths))./units_perDepth;
nonSomaUnits_perDepth = histcounts(template_depths(unitType == 3), min(template_depths):binSize:max(template_depths))./units_perDepth;
catch
    units_perDepth = histcounts(template_depths, min(template_depths):binSize:max(template_depths));
noiseUnits_perDepth = histcounts(template_depths, min(template_depths):binSize:max(template_depths))./units_perDepth;
goodUnits_perDepth = histcounts(template_depths, min(template_depths):binSize:max(template_depths))./units_perDepth;
nonSomaUnits_perDepth = histcounts(template_depths, min(template_depths):binSize:max(template_depths))./units_perDepth;

end
unittype_ax = subplot('Position', [0.4, 0.1, 0.1, 0.8]); hold on;
plot(smoothdata(goodUnits_perDepth, 'gaussian', [2 2]), min(template_depths)+binSize/2:binSize:max(template_depths)-binSize/2);
plot(smoothdata(noiseUnits_perDepth, 'gaussian', [2 2]), min(template_depths)+binSize/2:binSize:max(template_depths)-binSize/2);
plot(smoothdata(nonSomaUnits_perDepth, 'gaussian', [2 2]), min(template_depths)+binSize/2:binSize:max(template_depths)-binSize/2);
legend({'good', 'noise', 'non-soma'})
set(unittype_ax, 'YDir', 'reverse' )

%caxis([0,max(mua_corr(mua_corr ~= 1))]); colormap(hot);
ylim([min_depths, max_depths]);
set(unittype_ax, 'YTick', []);
title('Unit type');
set(unittype_ax, 'FontSize', 12)
xlabel(unittype_ax, 'Fraction');

% Plot LFP median-subtracted correlation
if length(lfp) > 1
    lfp_moving_median = 10; % channels to take sliding median
    lfp_ax = subplot('Position', [0.5, 0.1, 0.3, 0.8]);
    imagesc([min_depths, max_depths], [min_depths, max_depths], ...
        corrcoef((movmedian(zscore(double(lfp(theseChannels, :)), [], 2), lfp_moving_median, 1) - ...
        nanmedian(zscore(double(lfp(theseChannels, :)), [], 2), 1))'));
    xlim([min_depths, max_depths]);
    ylim([min_depths, max_depths]);
    set(lfp_ax, 'YTick', []);
    title('LFP power');
    set(lfp_ax, 'FontSize', 12)
    caxis([-1, 1])
    xlabel(lfp_ax, 'Depth (\mum)');
    colormap(lfp_ax, brewermap([], '*RdBu'));
end

% Plot probe areas (interactive)
% (load the colormap - located in the repository, find by associated fcn)
allenCCF_path = fileparts(which('allenCCFbregma'));
cmap_filename = [allenCCF_path, filesep, 'allen_ccf_colormap_2017.mat'];
load(cmap_filename);

probe_areas_ax = subplot('Position', [0.8, 0.1, 0.05, 0.8]);

% Convert probe CCF coordinates to linear depth (*25 to convert to um)
% (use the dorsal-most coordinate as the reference point)
[~, dv_sort_idx] = sort(probe_ccf(use_probe).trajectory_coords(:, 2));

probe_trajectory_depths = ...
    pdist2(probe_ccf(use_probe).trajectory_coords, ...
    probe_ccf(use_probe).trajectory_coords((dv_sort_idx == 1), :))*25;

trajectory_area_boundary_idx_non_linear = ...
    [1; find(diff(double(probe_ccf(use_probe).trajectory_areas)) ~= 0) + 1];
if trajectory_area_boundary_idx_non_linear(end) ~= size(probe_ccf(use_probe).trajectory_areas, 1)
    trajectory_area_boundary_idx_non_linear = [trajectory_area_boundary_idx_non_linear; size(probe_ccf(use_probe).trajectory_areas, 1)];
end
%trajectory_area_boundary_idx_non_linear = [trajectory_area_boundary_idx_non_linear; trajectory_area_boundary_idx_non_linear(end)+1];

% recreate linearily spaced areas 
% - depths
areas_linear_depth = min(probe_trajectory_depths):max(probe_trajectory_depths);
areas_linear = nan(size(areas_linear_depth, 1), 1, 1);

for iArea = 1:size(trajectory_area_boundary_idx_non_linear, 1) - 1
    theseDepths = find(areas_linear_depth >= probe_trajectory_depths(trajectory_area_boundary_idx_non_linear(iArea)) & ...
        areas_linear_depth <= probe_trajectory_depths(trajectory_area_boundary_idx_non_linear(iArea+1)));
    areas_linear(theseDepths) = probe_ccf(use_probe).trajectory_areas(trajectory_area_boundary_idx_non_linear(iArea)+1);
end
areas_linear = areas_linear';

% - boundaries
trajectory_area_boundary_idx = ...
    [1; find(diff(double(areas_linear)) ~= 0) + 1; size(areas_linear, 1)];

trajectory_area_boundaries = areas_linear_depth(trajectory_area_boundary_idx);

% - centers 
trajectory_area_centers = nan(size(trajectory_area_boundaries, 2)-1, 1);
for iBoundary = 1:size(trajectory_area_boundaries, 2) - 1
    trajectory_area_centers(iBoundary) = (trajectory_area_boundaries(iBoundary) + ...
        trajectory_area_boundaries(iBoundary+1)) ./ 2;
end

% -labels
trajectory_area_labels = cell(size(trajectory_area_boundaries, 2)-1, 1);
for iArea = 1:size(trajectory_area_boundaries, 2)
    trajectory_area_labels(iArea) = st.acronym(st.id == ...
        areas_linear(trajectory_area_boundary_idx(iArea)));
end

% plot
[~, area_dv_sort_idx] = sort(trajectory_area_centers);

probe_image = image([], areas_linear_depth, areas_linear);
colormap(probe_areas_ax, cmap);
caxis([1, size(cmap, 1)])
set(probe_areas_ax, 'YTick', (trajectory_area_centers(area_dv_sort_idx)-0.5), ...
    'YTickLabels', trajectory_area_labels(area_dv_sort_idx));
set(probe_areas_ax, 'XTick', []);
set(probe_areas_ax, 'YAxisLocation', 'right')

ylim([min_depths, max_depths]);
ylabel({'Probe areas', '(Arrow/shift keys to move)', '(Escape: save & quit)'});
set(probe_areas_ax, 'FontSize', 10)

% Draw boundary lines and make them draggable
boundary_lines = gobjects(length(trajectory_area_boundaries), 1);
for curr_boundary = 1:length(trajectory_area_boundaries)
    boundary_lines(curr_boundary) = line(probe_areas_ax, [-13.5, 1.5], ...
        repmat((trajectory_area_boundaries(curr_boundary)-0.5), 1, 2), 'color', 'b', ...
        'linewidth', 2, 'Tag', num2str(curr_boundary));
    draggable_line(boundary_lines(curr_boundary), probe_areas_ax, gui_fig);
end
set(probe_areas_ax, 'Clipping', 'off');

probe_ccf(use_probe).probe_depths_linear = areas_linear_depth;
probe_ccf(use_probe).probe_areas_linear = areas_linear;

% get unique region chunks and their scaling 
unique_cols = probe_image.CData(round(trajectory_area_centers(area_dv_sort_idx)));
scaling_per_region = ones(length(unique_cols),1);

% probe_depths (classic) 
gui_data.probe_trajectory_depths_simple = probe_trajectory_depths; 

% Package into gui
gui_data.probe_ccf_fn = probe_ccf_fn;
gui_data.probe_ccf = probe_ccf;
gui_data.use_probe = use_probe;
gui_data.unique_cols = unique_cols;
gui_data.scaling_per_region = scaling_per_region;

% axes
gui_data.unit_ax = unit_ax;
gui_data.multiunit_ax = multiunit_ax;
gui_data.probe_areas_ax = probe_areas_ax;
gui_data.probe_areas_ax_ylim = ylim(probe_areas_ax);
gui_data.probe_image = probe_image;

% trajectory and boundaries
gui_data.probe_trajectory_depths_non_linear = probe_trajectory_depths;
gui_data.probe_trajectory_depths = areas_linear_depth;
gui_data.original_region_boundaries = trajectory_area_boundaries;
gui_data.new_region_boundaries = gui_data.original_region_boundaries;
gui_data.areas_center = trajectory_area_centers(area_dv_sort_idx);
gui_data.min_depth = min_depths;
gui_data.max_depth = max_depths;
gui_data.trajectory_area_labels = trajectory_area_labels(area_dv_sort_idx);

% Upload gui data
guidata(gui_fig, gui_data);

end


function keypress(gui_fig, eventdata)

% Get guidata
gui_data = guidata(gui_fig);

% Set amounts to move by with/without shift
if any(strcmp(eventdata.Modifier, 'shift'))
    y_change = 100;
else
    y_change = 1;
end

switch eventdata.Key

    % up/down: move probe areas
    case 'uparrow'
        new_ylim = gui_data.probe_areas_ax_ylim - y_change;
        ylim(gui_data.probe_areas_ax, new_ylim);
        gui_data.probe_areas_ax_ylim = new_ylim;
        % Upload gui data
        guidata(gui_fig, gui_data);

    case 'downarrow'
        new_ylim = gui_data.probe_areas_ax_ylim + y_change;
        ylim(gui_data.probe_areas_ax, new_ylim);
        gui_data.probe_areas_ax_ylim = new_ylim;
        % Upload gui data
        guidata(gui_fig, gui_data);

       
    case 'escape' % escape: save and quit
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?', 'Confirm exit', opts);
        if strcmp(user_confirm, 'Yes')

            probe_ccf = gui_data.probe_ccf;

            % shift the non-linear to the linear 
            iProbe = gui_data.use_probe;
          
            probe_ccf(iProbe).probe_depths = gui_data.probe_trajectory_depths_simple;
            probe_ccf(iProbe).probe_depths =  probe_ccf(iProbe).probe_depths - probe_ccf(iProbe).probe_depths(1) ...
                 - gui_data.probe_areas_ax_ylim(1);

            % area boundary depths
            area_boundary_depth = gui_data.original_region_boundaries - gui_data.probe_areas_ax_ylim(1);

            % SCALING FACTORS 
            ori_sizes = diff(gui_data.original_region_boundaries);
            new_sizes = diff(gui_data.new_region_boundaries);
            scaling = new_sizes./ori_sizes;
            % vector of depths (linear)
            depth_diff_vector = diff(probe_ccf(iProbe).probe_depths);

            % vector of scale factors 
            depth_scale_vector = ones(size(probe_ccf(iProbe).probe_depths,1),1);
            for iArea = 1:size(gui_data.scaling_per_region, 1)
                start=find(probe_ccf(iProbe).probe_depths >area_boundary_depth(iArea), 1,'first');
                stop=find(probe_ccf(iProbe).probe_depths >area_boundary_depth(iArea+1), 1,'first');

                depth_scale_vector(start:stop)...
                    = scaling(iArea);
            end

            % new depths: value 1 + (diff * scale vectors)
            scaled_depths = [probe_ccf(iProbe).probe_depths(1);...
                probe_ccf(iProbe).probe_depths(1) + cumsum(depth_diff_vector .* depth_scale_vector(2:end))];
            
            % % sanity checks 
            % probe_ccf(iProbe).probe_depths(find(probe_ccf(iProbe).trajectory_areas==1022, 1, 'first'))
            % 
            % probe_ccf(iProbe).probe_depths(find(probe_ccf(iProbe).trajectory_areas==6, 1, 'first'))
            % probe_ccf(iProbe).probe_depths(find(probe_ccf(iProbe).trajectory_areas==6, 1, 'last'))
            % 
            % probe_ccf(iProbe).probe_depths(find(probe_ccf(iProbe).trajectory_areas==672, 1, 'first'))
            % probe_ccf(iProbe).probe_depths(find(probe_ccf(iProbe).trajectory_areas==672, 1, 'last'))
            % 
            %  % sanity checks 
            % scaled_depths(find(probe_ccf(iProbe).trajectory_areas==1022, 1, 'first'))
            % 
            % scaled_depths(find(probe_ccf(iProbe).trajectory_areas==6, 1, 'first'))
            % scaled_depths(find(probe_ccf(iProbe).trajectory_areas==6, 1, 'last'))
            % 
            % scaled_depths(find(probe_ccf(iProbe).trajectory_areas==672, 1, 'first'))
            % scaled_depths(find(probe_ccf(iProbe).trajectory_areas==672, 1, 'last'))

            probe_ccf(iProbe).probe_depths=scaled_depths;


            % previous: probe_ccf(iProbe).probe_depths =  probe_ccf(iProbe).probe_depths - gui_data.probe_areas_ax_ylim(1);

            
            % region boundaries 
            probe_ccf(iProbe).trajectory_coords_aligned = probe_ccf(iProbe).trajectory_coords;
            probe_ccf(iProbe).trajectory_areas_aligned = probe_ccf(iProbe).trajectory_areas;
           
            
            % save scaling
            probe_ccf(iProbe).area_scaling = gui_data.scaling_per_region;

            % Save the appended probe_ccf structure
            save(gui_data.probe_ccf_fn, 'probe_ccf');

            % Close the figure
            close(gui_fig);
        end


end

end


function draggable_line(hLine, ax, fig)
set(hLine, 'ButtonDownFcn', @startDragFcn);
% Get guidata
gui_data = guidata(fig);


    function startDragFcn(src, ~)
        src.UserData = get(ax, 'CurrentPoint');
        set(fig, 'WindowButtonMotionFcn', {@draggingFcn, src, ax, fig});
        set(fig, 'WindowButtonUpFcn', @stopDragFcn);
    end

    function draggingFcn(~, ~, src, ax, fig)
        curr_pt = get(ax, 'CurrentPoint');
        orig_pt = src.UserData;

        dy = curr_pt(1, 2) - orig_pt(1, 2);
        ydata = get(src, 'YData');
        ydata = ydata + dy;
        set(src, 'YData', ydata);

        src.UserData = curr_pt;
        update_probe_areas_image(ax, src, fig);
    end

    function stopDragFcn(~, ~)
        set(fig, 'WindowButtonMotionFcn', '');
        set(fig, 'WindowButtonUpFcn', '');

    end

    function update_probe_areas_image(probe_areas_ax, src, fig)
        % (get info)
        newLine_location = src.UserData(1, 2);
        lineTag = str2num(get(src, 'Tag'));
        gui_data = guidata(fig);
       % original_depth = gui_data.probe_trajectory_depths(gui_data.original_region_boundaries(lineTag));
        prev_depth = gui_data.probe_trajectory_depths(gui_data.new_region_boundaries(lineTag)); % up or down?
        depth_change_stretch = newLine_location > prev_depth;
        
        % (get new area boundaries)
        new_colors = gui_data.probe_image.CData;
        gui_data.new_region_boundaries = gui_data.original_region_boundaries;
        gui_data.new_region_boundaries(lineTag) = find(gui_data.probe_trajectory_depths >= newLine_location, 1, 'first');

        % update scaling factor 
        original_size_up = gui_data.probe_trajectory_depths(gui_data.original_region_boundaries(lineTag)) -...
            gui_data.probe_trajectory_depths(gui_data.original_region_boundaries(lineTag-1));
        new_size_up = gui_data.probe_trajectory_depths(gui_data.new_region_boundaries(lineTag)) -...
            gui_data.probe_trajectory_depths(gui_data.new_region_boundaries(lineTag-1));
        gui_data.scaling_per_region(lineTag-1) = new_size_up/original_size_up;
        
        original_size_down = gui_data.probe_trajectory_depths(gui_data.original_region_boundaries(lineTag+1)) -...
            gui_data.probe_trajectory_depths(gui_data.original_region_boundaries(lineTag));
        new_size_down = gui_data.probe_trajectory_depths(gui_data.new_region_boundaries(lineTag+1)) -...
            gui_data.probe_trajectory_depths(gui_data.new_region_boundaries(lineTag));
        gui_data.scaling_per_region(lineTag) = new_size_down/original_size_down;
    
        % (update image)
        if depth_change_stretch == 1 % stretch up
            curr_area_color = gui_data.unique_cols(lineTag-1);
            new_colors(round(gui_data.new_region_boundaries(lineTag-1)+1):find(gui_data.probe_trajectory_depths >= newLine_location, 1, 'first'))...
                = curr_area_color;
        else % stretch down
            curr_area_color = gui_data.unique_cols(lineTag);
            new_colors(find(gui_data.probe_trajectory_depths >= newLine_location, 1, 'first'):round(gui_data.areas_center(lineTag)))...
                = curr_area_color;
        end
        set(gui_data.probe_image, 'CData', new_colors) % update color

        % (update centers)
        gui_data.area_centers = nan(size(gui_data.new_region_boundaries, 2)-1, 1);
        for iBoundary = 1:size(gui_data.new_region_boundaries, 2) - 1
            gui_data.area_centers(iBoundary) = (gui_data.new_region_boundaries(iBoundary) + ...
                gui_data.new_region_boundaries(iBoundary+1)) ./ 2;
        end
        set(probe_areas_ax, 'YTick', (gui_data.area_centers -0.5), ...
            'YTickLabels', gui_data.trajectory_area_labels);

        % (update coordinates) + ineed to update this when shifting, no
        % longer only in saving. or just have an updates stretch factpor
        % thing????
        %gui_data.


        % Upload gui data
        guidata(fig, gui_data);

    end
end
