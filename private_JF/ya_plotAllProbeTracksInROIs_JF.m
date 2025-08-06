function ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsOfInterest, patchBrain, onlyROIProbes, showPoints, useBezierFit, showRegionPlot, regionColors, blackBackground, thickBrainLines)
% ya_plotAllProbeTracksInROIs_JF - Plot probe tracks for multiple animals with regions of interest
%
% Inputs:
%   theseAnimals - Cell array of animal names (e.g., {'JF058', 'JF059'})
%   regionsOfInterest - Cell array of region names to highlight (e.g., {'CP', 'GPe', 'SNr'})
%                       If empty or not provided, defaults to {'CP', 'GPe', 'GPi', 'STN', 'SNr'}
%   patchBrain - Boolean to use surface patch (1) or wire grid (0) for brain. Default: 0
%   onlyROIProbes - Boolean to plot only probes that pass through ROIs (1) or all probes (0). Default: 0
%   showPoints - Boolean to show probe points (1) or just fitted lines (0). Default: 1
%   useBezierFit - Boolean to use Bezier curve fit (1) or linear fit (0). Default: 1
%   showRegionPlot - Boolean to show brain regions per probe plot (1) or not (0). Default: 1
%   regionColors - Cell array of RGB colors for each region, or 'allen' to use Allen CCF colormap. Default: 'allen'
%   blackBackground - Boolean to use black background (1) or white (0). Default: 0
%   thickBrainLines - Line width for brain grid (default: 0.5, thicker: 2.0). Default: 0.5
%
% Example:
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 1, 1, 1, 1); % All features with Allen colors
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 0, 0, 1, 0, {[1 0 0], [0 1 0]}); % Custom colors
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 0, 0, 1, 1, 'allen', 1, 2); % Black background, thick lines

% Set defaults
if nargin < 2 || isempty(regionsOfInterest)
    regionsOfInterest = {'CP', 'GPe', 'GPi', 'STN', 'SNr'};
end
if nargin < 3
    patchBrain = 0;
end
if nargin < 4
    onlyROIProbes = 0;
end
if nargin < 5
    showPoints = 1;
end
if nargin < 6
    useBezierFit = 1;
end
if nargin < 7
    showRegionPlot = 1;
end
if nargin < 8 || isempty(regionColors)
    regionColors = 'allen';
end
if nargin < 9
    blackBackground = 0;
end
if nargin < 10
    thickBrainLines = 0.5;
end

% Initialize paths and parameters
cl_myPaths;
% Bregma will be set after loading atlas
animalsType = {'Naive'};
regionsNames = regionsOfInterest;
regions = regionsOfInterest;
regionPlotLoc = repmat([-1], 1, length(regionsOfInterest)); % Default to left hemisphere

% Load Allen atlas - use 10um version as expected by the original code
allen_atlas_path = '/home/jf5479/Dropbox/Atlas/allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
st = ya_loadStructureTree([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
slice_spacing = 10;

% Also load the 25um atlas for st_br compatibility
brainglobeLocation = '/home/jf5479/Dropbox/Atlas/brainglobe/';
atlasType = 'allen';
atlasSpecies = 'mouse';
atlasResolution_um = 25;
atlasLocation = dir([brainglobeLocation, atlasType, '_', atlasSpecies, '_', num2str(atlasResolution_um), 'um*']);
[~, ~, st_br, ~] = ya_loadAllenAtlas([atlasLocation.folder, filesep, atlasLocation.name]);
n_str = 0;
n_gpe =0;
n_snr =0;

% Generate colors for regions
if ischar(regionColors) && strcmp(regionColors, 'allen')
    % Load Allen CCF colormap
    cmap_filename = [allen_atlas_path, filesep, 'allen_ccf_colormap_2017.mat'];
    if exist(cmap_filename, 'file')
        load(cmap_filename, 'cmap');
        theseColors = cell(length(regionsNames), 1);
        
        % Need to use st_br for proper ID mapping as in trajectory_areas
        for iRegion = 1:length(regionsNames)
            % Find structure in st_br (brainglobe structure tree) which matches trajectory_areas
            struct_idx_br = find(strcmp(st_br.acronym, regionsNames{iRegion}));
            if ~isempty(struct_idx_br)
                % Get the structure ID from st_br - this matches trajectory_areas values
                structure_id = st_br.id(struct_idx_br(1));
                
                % The colormap is indexed by structure ID directly
                if structure_id > 0 && structure_id <= size(cmap, 1)
                    theseColors{iRegion} = cmap(structure_id, :);
                else
                    % Fallback color
                    theseColors{iRegion} = [0.5, 0.5, 0.5];
                end
            else
                % If not found in st_br, try st (Allen structure tree)
                struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
                if ~isempty(struct_idx)
                    structure_id = st.id(struct_idx(1));
                    if structure_id > 0 && structure_id <= size(cmap, 1)
                        theseColors{iRegion} = cmap(structure_id, :);
                    else
                        theseColors{iRegion} = [0.5, 0.5, 0.5];
                    end
                else
                    theseColors{iRegion} = [0.5, 0.5, 0.5];
                end
            end
        end
    else
        % Fallback if colormap not found
        warning('Allen CCF colormap not found, using default colors');
        theseColors = lines(length(regionsNames));
        theseColors = mat2cell(theseColors, ones(size(theseColors,1),1), 3);
    end
elseif iscell(regionColors) && length(regionColors) == length(regionsNames)
    % User-provided colors
    theseColors = regionColors;
else
    % Default to lines colormap
    theseColors = lines(length(regionsNames));
    theseColors = mat2cell(theseColors, ones(size(theseColors,1),1), 3);
end
%add probe types, depths
for iType = 1:size(animalsType, 2)
    %theseTypes = strcmp(recordingInfo.Type, animalsType{iType});
    %theseColors = {rgb('DeepSkyBlue'); rgb('DarkOrange'); rgb('Hotpink'); rgb('SeaGreen'); rgb('Crimson')};

    slice_spacing = 10;
    structure_alpha = 0.2;
    %get colors (overrride allen)
    %figure();
    if patchBrain
        ya_plotBrainSurface(allenAtlasPath)
    else
        if blackBackground
            [~, brain_outline] = plotBrainGrid([], [], [], true);
        else
            [~, brain_outline] = plotBrainGrid([], []);
        end
        % Set line width for brain grid
        set(brain_outline, 'LineWidth', thickBrainLines);
    end

    %overlay regions - plot bilaterally (both hemispheres)
    for iRegion = 1:length(regionsNames)
        curr_plot_structure = find(strcmp(st.acronym, regionsNames{iRegion}));
        
        % Plot the full structure (both hemispheres)
        structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);
        
        if strcmp(regionsNames{iRegion}, 'STR') == 0 && ~isempty(structure_3d.vertices)
            hold on;
            axis vis3d equal off manual
            view([-30, 25]);
            caxis([0, 300]);
            [ap_max, dv_max, ml_max] = size(tv);
            xlim([-10, ap_max + 10])
            ylim([-10, ml_max + 10])
            zlim([-10, dv_max + 10])
            structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
                'Faces', structure_3d.faces, ...
                'FaceColor', theseColors{iRegion, :}, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
        end
    end % end region plotting loop
    
    % Plot probe tracks for all animals
    % Get colors for each animal (same color per mouse)
    mouseColors = lines(length(theseAnimals));
    
    for iAnimal = 1:size(theseAnimals, 2)
            %iAnimal = iAnimal + 1;
            % Load probe data using the same path pattern as ya_histologyMain_JF_wittenLab
            outputDir = ['/home/jf5479/cup/Chris/data/cta_backwards/' theseAnimals{iAnimal} '/histology/alignedAllen/'];
            probe_ccf_location = [outputDir, 'probe_ccf.mat'];
            
            % Check if file exists
            if ~exist(probe_ccf_location, 'file')
                warning('probe_ccf.mat not found for %s at %s', theseAnimals{iAnimal}, probe_ccf_location);
                continue;
            end
            
            load(probe_ccf_location)

            % Get color for this animal (same for all probes from this mouse)
            animalColor = mouseColors(iAnimal, :);
            
            for iProbe = 1:size(probe_ccf, 1)
                curr_probe = iProbe;
                
                % Check if probe has points
                if isfield(probe_ccf(iProbe), 'points') && ~isempty(probe_ccf(iProbe).points)
                    
                    % Check if we should plot this probe based on ROI filtering
                    plotThisProbe = true;
                    if onlyROIProbes
                        plotThisProbe = false;
                        % Only plot if probe passes through any of the ROIs
                        if isfield(probe_ccf(iProbe), 'trajectory_areas')
                            for iRegion = 1:length(regionsNames)
                                struct_curr_br = st_br.id(strcmp(st_br.acronym, regionsNames{iRegion}));
                                if any(probe_ccf(iProbe).trajectory_areas == struct_curr_br)
                                    plotThisProbe = true;
                                    break;
                                end
                            end
                        end
                    end
                    
                    if plotThisProbe
                        thesePoints = probe_ccf(curr_probe).points * 2.5; % Scale to match brain grid
                        
                        % Plot probe points if requested
                        if showPoints
                            plot3(thesePoints(:, 1), ...
                                thesePoints(:, 2), ...
                                thesePoints(:, 3), ...
                                '.', 'color', animalColor, 'MarkerSize', 20);
                        end
                        
                        % Fit curve through points
                        if useBezierFit && size(thesePoints, 1) >= 3
                            % Use Bezier curve fitting
                            % Sort points by one dimension (e.g., DV) to get proper curve order
                            [~, sort_idx] = sort(thesePoints(:, 2)); % Sort by DV (y-axis)
                            sorted_points = thesePoints(sort_idx, :);
                            
                            % Create Bezier curve
                            t = linspace(0, 1, 1000);
                            bezier_curve_points = bezier_curve(t, sorted_points);
                            
                            % Plot Bezier curve
                            plot3(bezier_curve_points(:, 1), bezier_curve_points(:, 2), bezier_curve_points(:, 3), ...
                                'color', animalColor, 'linewidth', 2);
                        else
                            % Use linear fit (original method)
                            r0 = mean(thesePoints, 1);
                            xyz = bsxfun(@minus, thesePoints, r0);
                            [~, ~, V] = svd(xyz, 0);
                            histology_probe_direction = V(:, 1);
                            
                            % Make sure the direction goes down in DV - flip if it's going up
                            if histology_probe_direction(2) < 0
                                histology_probe_direction = -histology_probe_direction;
                            end

                            line_eval = [-1000, 1000];
                            probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0);
                            
                            % Plot linear fit
                            plot3(probe_fit_line(:, 1), probe_fit_line(:, 2), probe_fit_line(:, 3), ...
                                'color', animalColor, 'linewidth', 2);
                        end
                    end
                    
                    % Count probes per region (regardless of plotting)
                    if isfield(probe_ccf(iProbe), 'trajectory_areas')
                        for iRegion = 1:length(regionsNames)
                            struct_curr_br = st_br.id(strcmp(st_br.acronym, regionsNames{iRegion}));
                            if any(probe_ccf(iProbe).trajectory_areas == struct_curr_br)
                                % Count probes per region
                                if strcmpi(regionsNames{iRegion}, 'CP') || strcmpi(regionsNames{iRegion}, 'STR')
                                    n_str = n_str+1;
                                elseif strcmpi(regionsNames{iRegion}, 'GPe')
                                    n_gpe = n_gpe+1;
                                elseif strcmpi(regionsNames{iRegion}, 'SNr')
                                    n_snr = n_snr+1;
                                end
                            end
                        end
                    end
                end
            end
    end % end animal loop
end % end animalsType loop
% Set final view
view([-30, 25]);
if blackBackground
    set(gcf, 'Color', 'k');
    set(gca, 'Color', 'k');
    % Update axes colors for black background
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.ZColor = 'w';
else
    set(gcf, 'Color', 'w');
end

% Enable 3D rotation
h = rotate3d(gca);
h.Enable = 'on';

% Add title
t = title(sprintf('Probe tracks for %d animals in %s', length(theseAnimals), strjoin(regionsNames, ', ')));
if blackBackground
    set(t, 'Color', 'w');
end

% Add legends for both mouse colors and region colors
legendEntries = {};
legendHandles = [];

% Mouse legend entries
for iMouse = 1:length(theseAnimals)
    % Create a dummy line for legend
    h_legend = plot3(NaN, NaN, NaN, 'color', mouseColors(iMouse, :), 'linewidth', 3);
    legendHandles(end+1) = h_legend;
    legendEntries{end+1} = strrep(theseAnimals{iMouse},'_', '-');
end

% Add separator
h_legend = plot3(NaN, NaN, NaN, 'color', 'none');
legendHandles(end+1) = h_legend;
legendEntries{end+1} = '--- Regions ---';

% Region legend entries
for iRegion = 1:length(regionsNames)
    % Create a patch for region legend
    h_patch = patch(NaN, NaN, NaN, theseColors{iRegion}, 'EdgeColor', 'none');
    legendHandles(end+1) = h_patch;
    legendEntries{end+1} = regionsNames{iRegion};
end

% Create legend
if ~isempty(legendHandles)
    leg = legend(legendHandles, legendEntries, 'Location', 'bestoutside', 'FontSize', 10);
    if blackBackground
        set(leg, 'TextColor', 'w', 'Color', 'k', 'EdgeColor', 'w');
    end
end

%% Create brain regions per probe plot
if showRegionPlot
    plotProbeRegions(theseAnimals, mouseColors, onlyROIProbes, regionsNames, st_br, allen_atlas_path);
end

end

function B = bezier_curve(t, control_points)
% Bezier curve function
n = size(control_points, 1) - 1; % degree of the polynomial
B = zeros(3, length(t)); % 3D curve

[~, control_points_idx] = sort(control_points(:, 3));
control_points = control_points(control_points_idx, :);
for i = 0:n
    B = B + nchoosek(n, i) * (1 - t).^(n - i) .* t.^i .* control_points(i+1, :)';
end
B = B';
end

function plotProbeRegions(theseAnimals, mouseColors, onlyROIProbes, regionsNames, st_br, allen_atlas_path)
% Plot brain regions that each probe passes through
% Load Allen CCF colormap
cmap_filename = [allen_atlas_path, filesep, 'allen_ccf_colormap_2017.mat'];
if ~exist(cmap_filename, 'file')
    warning('Allen CCF colormap not found, using default colormap');
    cmap = lines(256);
else
    load(cmap_filename, 'cmap');
end

% Count total probes across all animals
totalProbes = 0;
allProbeData = {};

for iAnimal = 1:length(theseAnimals)
    % Load probe data
    outputDir = ['/home/jf5479/cup/Chris/data/cta_backwards/' theseAnimals{iAnimal} '/histology/alignedAllen/'];
    probe_ccf_location = [outputDir, 'probe_ccf.mat'];
    
    if ~exist(probe_ccf_location, 'file')
        continue;
    end
    
    load(probe_ccf_location, 'probe_ccf');
    
    for iProbe = 1:length(probe_ccf)
        % Check if we should include this probe based on ROI filtering
        includeThisProbe = true;
        if onlyROIProbes
            includeThisProbe = false;
            if isfield(probe_ccf(iProbe), 'trajectory_areas') && ~isempty(probe_ccf(iProbe).trajectory_areas)
                for iRegion = 1:length(regionsNames)
                    struct_curr_br = st_br.id(strcmp(st_br.acronym, regionsNames{iRegion}));
                    if any(probe_ccf(iProbe).trajectory_areas == struct_curr_br)
                        includeThisProbe = true;
                        break;
                    end
                end
            end
        end
        
        if includeThisProbe && isfield(probe_ccf(iProbe), 'trajectory_areas') && ~isempty(probe_ccf(iProbe).trajectory_areas)
            totalProbes = totalProbes + 1;
            allProbeData{totalProbes} = struct(...
                'trajectory_areas', probe_ccf(iProbe).trajectory_areas, ...
                'animal', theseAnimals{iAnimal}, ...
                'probe_num', iProbe, ...
                'animal_idx', iAnimal);
        end
    end
end

if totalProbes == 0
    warning('No probes with trajectory data found');
    return;
end

% Create new figure for region plots
figure('Name', 'Brain Regions per Probe', 'Color', 'k', 'Position', [100, 100, max(1200, totalProbes*100), 600]);

% Create subplots
for iProbe = 1:totalProbes
    curr_axes = subplot(1, totalProbes, iProbe);
    
    trajectory_areas = allProbeData{iProbe}.trajectory_areas;
    animal_name = allProbeData{iProbe}.animal;
    probe_num = allProbeData{iProbe}.probe_num;
    animal_idx = allProbeData{iProbe}.animal_idx;
    
    % Find area boundaries and centers
    trajectory_area_boundaries = [1; find(diff(trajectory_areas) ~= 0); length(trajectory_areas)];
    trajectory_area_centers = trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries) / 2;
    
    if length(trajectory_area_centers) > 1
        % Get area labels
        trajectory_area_labels = cell(length(trajectory_area_centers), 1);
        for iArea = 1:length(trajectory_area_centers)
            area_id = trajectory_areas(round(trajectory_area_centers(iArea)));
            matching_idx = find(st_br.id == area_id);
            if ~isempty(matching_idx)
                trajectory_area_labels{iArea} = st_br.acronym{matching_idx(1)};
            else
                trajectory_area_labels{iArea} = sprintf('ID_%d', area_id);
            end
        end
        
        % Plot trajectory areas
        image(trajectory_areas);
        colormap(curr_axes, cmap);
        caxis([1, size(cmap, 1)]);
        set(curr_axes, 'YTick', trajectory_area_centers, 'YTickLabels', trajectory_area_labels);
        set(curr_axes, 'XTick', []);
        
        % Style the plot
        set(curr_axes, 'Color', 'k');
        set(curr_axes, 'YColor', 'w');
        set(curr_axes, 'XColor', 'w');
        set(curr_axes, 'GridColor', 'w');
        
        % Set title with animal and probe info, colored by mouse
        title_str = sprintf('%s-P%d', strrep(animal_name, '_', '-'), probe_num);
        title(title_str, 'Color', mouseColors(animal_idx, :), 'FontSize', 10);
    else
        % Empty plot for probes without trajectory data
        set(curr_axes, 'Color', 'k');
        axis off;
        title(sprintf('%s-P%d\n(No data)', strrep(animal_name, '_', '-'), probe_num), ...
            'Color', mouseColors(animal_idx, :), 'FontSize', 10);
    end
end

% Add overall title
sgtitle('Brain Regions Traversed by Each Probe', 'FontSize', 14, 'Color', 'w');

end

