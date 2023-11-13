function ya_alignEphysAndHistology(st, outputDir, ...
    spike_times, spike_templates, template_depths, spike_xdepths, template_xdepths, ...
    lfp, lfp_channel_positions, lfp_channel_xpositions, use_probe, isSpikeGlx, curr_shank)
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
if isSpikeGlx
    n_corr_groups = 20;
else
    n_corr_groups = 40;
end
if isSpikeGlx && max(template_depths) <= 384 * 7.5 %2.0 probes, shorter shank
    max_depths = 384 * 7.5;
    min_depths = 0;
    depth_group_edges = linspace(min_depths, max_depths, n_corr_groups+1);
else
    min_depths = 0;
    max_depths = 3840; % (hardcode, sometimes kilosort2 drops channels)
    depth_group_edges = linspace(0, max_depths, n_corr_groups+1);
end

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
multiunit_ax = subplot('Position', [0.2, 0.1, 0.3, 0.8]);
imagesc(depth_group_centers, depth_group_centers, mua_corr);
%caxis([0,max(mua_corr(mua_corr ~= 1))]); colormap(hot);
ylim([min_depths, max_depths]);
set(multiunit_ax, 'YTick', []);
title('MUA correlation');
set(multiunit_ax, 'FontSize', 12)
xlabel(multiunit_ax, 'Multiunit depth');

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
    probe_ccf(use_probe).trajectory_coords((dv_sort_idx == 1), :));

%probe_trajectory_depths = probe_ccf(use_probe).trajectory_coords(:,2);

trajectory_area_boundary_idx_non_linear = ...
    [1; find(diff(double(probe_ccf(use_probe).trajectory_areas)) ~= 0) + 1];
if trajectory_area_boundary_idx_non_linear(end) ~= size(probe_ccf(use_probe).trajectory_areas, 1)
    trajectory_area_boundary_idx_non_linear = [trajectory_area_boundary_idx_non_linear; size(probe_ccf(use_probe).trajectory_areas, 1)];
end
% recreate linearily spaced areas
areas_linear_depth = min(probe_trajectory_depths):max(probe_trajectory_depths);

areas_linear = nan(size(areas_linear_depth, 1), 1, 1);

for iBoundary = 1:size(trajectory_area_boundary_idx_non_linear, 1) - 1

    theseDepths = find(areas_linear_depth >= probe_trajectory_depths(trajectory_area_boundary_idx_non_linear(iBoundary)) & ...
        areas_linear_depth <= probe_trajectory_depths(trajectory_area_boundary_idx_non_linear(iBoundary+1)));

    areas_linear(theseDepths) = probe_ccf(use_probe).trajectory_areas(trajectory_area_boundary_idx_non_linear(iBoundary)+1);

end
areas_linear = areas_linear';
trajectory_area_boundary_idx = ...
    [1; find(diff(double(areas_linear)) ~= 0) + 1; size(areas_linear, 1)];


trajectory_area_boundaries = areas_linear_depth(trajectory_area_boundary_idx);

trajectory_area_centers = nan(size(trajectory_area_boundaries, 2)-1, 1);
for iBoundary = 1:size(trajectory_area_boundaries, 2) - 1
    trajectory_area_centers(iBoundary) = (trajectory_area_boundaries(iBoundary) + ...
        trajectory_area_boundaries(iBoundary+1)) ./ 2;
end


trajectory_area_labels = cell(size(trajectory_area_boundaries, 2)-1, 1);
for iArea = 1:size(trajectory_area_boundaries, 2)
    trajectory_area_labels(iArea) = st.acronym(st.id == ...
        areas_linear(trajectory_area_boundary_idx(iArea)));
end
%trajectory_area_labels = st.acronym(st.id ==areas_linear(trajectory_area_boundary_idx));

[~, area_dv_sort_idx] = sort(trajectory_area_centers);

image([], areas_linear_depth*25, areas_linear);
colormap(probe_areas_ax, cmap);
caxis([1, size(cmap, 1)])
set(probe_areas_ax, 'YTick', trajectory_area_centers(area_dv_sort_idx)*25, ...
    'YTickLabels', trajectory_area_labels(area_dv_sort_idx));
set(probe_areas_ax, 'XTick', []);
set(probe_areas_ax, 'YAxisLocation', 'right')

ylim([min_depths, max_depths]);
ylabel({'Probe areas', '(Arrow/shift keys to move)', '(Escape: save & quit)'});
set(probe_areas_ax, 'FontSize', 10)

% Draw boundary lines at borders (and undo clipping to extend across all)
boundary_lines = gobjects;
for curr_boundary = 1:length(trajectory_area_boundaries)
    boundary_lines(curr_boundary, 1) = line(probe_areas_ax, [-13.5, 1.5], ...
        repmat((trajectory_area_boundaries(curr_boundary) - 0.5)*25, 1, 2), 'color', 'b', 'linewidth', 2);
end
set(probe_areas_ax, 'Clipping', 'off');

probe_ccf(use_probe).probe_depths_linear = areas_linear_depth;
probe_ccf(use_probe).probe_areas_linear = areas_linear;

% Package into gui
gui_data = struct;
gui_data.probe_ccf_fn = probe_ccf_fn;

gui_data.probe_ccf = probe_ccf;
gui_data.use_probe = use_probe;

gui_data.unit_ax = unit_ax;
gui_data.multiunit_ax = multiunit_ax;

gui_data.probe_areas_ax = probe_areas_ax;
gui_data.probe_areas_ax_ylim = ylim(probe_areas_ax);
gui_data.probe_trajectory_depths_non_linear = probe_trajectory_depths;

gui_data.probe_trajectory_depths = areas_linear_depth*25;

% Upload gui data
guidata(gui_fig, gui_data);

end


function keypress(gui_fig, eventdata)

% Get guidata
gui_data = guidata(gui_fig);

% Set amounts to move by with/without shift
if any(strcmp(eventdata.Modifier, 'shift'))
    y_change = 100;
    s_change = 0.1;
else
    y_change = 1;
    s_change = 0.01;
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

        % escape: save and quit
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?', 'Confirm exit', opts);
        if strcmp(user_confirm, 'Yes')

            probe_ccf = gui_data.probe_ccf;

            % Get the probe depths corresponding to the trajectory areas
            probe_depths_linear = gui_data.probe_trajectory_depths - ...
                gui_data.probe_areas_ax_ylim(1);

            probe_ccf(gui_data.use_probe).probe_depths_linear = probe_depths_linear;

            % Save the appended probe_ccf structure
            save(gui_data.probe_ccf_fn, 'probe_ccf');

            % Close the figure
            close(gui_fig);
        end
    case 'n' % stretch
        new_ylim = [gui_data.probe_areas_ax_ylim(1), gui_data.probe_areas_ax_ylim(2) .* (1 - s_change)];
        ylim(gui_data.probe_areas_ax, new_ylim);
        gui_data.probe_areas_ax_ylim = new_ylim;
        % Upload gui data
        guidata(gui_fig, gui_data);

    case 'w' % stretch
        probe_ccf = gui_data.probe_ccf;

        % Get the probe depths corresponding to the trajectory areas
        probe_depths_linear = gui_data.probe_trajectory_depths_linear .* (1 + s_change);

        probe_ccf(gui_data.use_probe).probe_depths_linear = probe_depths_linear;

        % Save the appended probe_ccf structure
        save(gui_data.probe_ccf_fn, 'probe_ccf');


end

end
