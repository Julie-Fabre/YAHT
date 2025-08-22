function ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsOfInterest, patchBrain, onlyROIProbes, showPoints, useBezierFit, showRegionPlot, regionColors, blackBackground, thickBrainLines, plotHemisphere)
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
%   plotHemisphere - Array specifying hemisphere for each region: -1 for left, 1 for right, 0 for both.
%                    If a probe passes through multiple regions with different hemispheres, it will be plotted multiple times.
%                    Default: empty (plots probes in their actual positions)
%
% Example:
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 1, 1, 1, 1); % All features with Allen colors
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 0, 0, 1, 0, {[1 0 0], [0 1 0]}); % Custom colors
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 0, 0, 1, 1, 'allen', 1, 2); % Black background, thick lines
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 0, 0, 1, 1, 'allen', 0, 0.5, [-1, 1]); % CP on left, SNr on right

% Set defaults
if nargin < 2 || isempty(regionsOfInterest)
    regionsOfInterest = {'CP', 'GPe', 'SNr'};
end
if nargin < 3 || isempty(patchBrain)
    patchBrain = 0;
end
if nargin < 4 || isempty(onlyROIProbes)
    onlyROIProbes = 0;
end
if nargin < 5 || isempty(showPoints)
    showPoints = 1;
end
if nargin < 6 || isempty(useBezierFit)
    useBezierFit = 1;
end
if nargin < 7 || isempty(showRegionPlot)
    showRegionPlot = 1;
end
if nargin < 8 || isempty(regionColors)
    regionColors = 'allen';
end
if nargin < 9 || isempty(blackBackground)
    blackBackground = 0;
end
if nargin < 10 || isempty(thickBrainLines)
    thickBrainLines = 0.5;
end
if nargin < 11 || isempty(plotHemisphere)
    plotHemisphere = [];
end

% Initialize paths and parameters
cl_myPaths;
% Bregma will be set after loading atlas
animalsType = {'Naive'};
regionsNames = regionsOfInterest;
regions = regionsOfInterest;
% Set up hemisphere plotting preferences
if ~isempty(plotHemisphere)
    if length(plotHemisphere) ~= length(regionsOfInterest)
        error('plotHemisphere must have the same length as regionsOfInterest');
    end
    regionPlotLoc = plotHemisphere;
else
    regionPlotLoc = zeros(1, length(regionsOfInterest)); % Default to both hemispheres (no mirroring)
end

% Load Allen atlas - use 10um version as expected by the original code
cl_myPaths;
%allenAtlasPath = '/home/jf5479/Dropbox/Atlas/allenCCF';
tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']);
av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']);
st = ya_loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']);
slice_spacing = 10;

% Also load the 25um atlas for st_br compatibility
%brainglobeLocation = '/home/jf5479/Dropbox/Atlas/brainglobe/';
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
    cmap_filename = [allenAtlasPath, filesep, 'allen_ccf_colormap_2017.mat'];
    if exist(cmap_filename, 'file')
        load(cmap_filename, 'cmap');
        theseColors = cell(length(regionsNames), 1);
        
        % Need to use Allen CCF (st) for proper ID mapping as trajectory_areas uses Allen IDs
        for iRegion = 1:length(regionsNames)
            % Find structure in Allen CCF structure tree which matches trajectory_areas
            struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
            if ~isempty(struct_idx)
                % Get the structure ID from Allen CCF - this matches trajectory_areas values
                structure_id = st.id(struct_idx(1));
                
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
    % User-provided colors as cell array
    theseColors = regionColors;
elseif isnumeric(regionColors) && size(regionColors, 1) == length(regionsNames) && size(regionColors, 2) == 3
    % User-provided colors as numeric matrix - convert to cell array
    theseColors = mat2cell(regionColors, ones(size(regionColors,1),1), 3);
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

    %overlay regions - plot based on hemisphere preferences
    for iRegion = 1:length(regionsNames)
        curr_plot_structure = find(strcmp(st.acronym, regionsNames{iRegion}));
        
        % Get the structure volume
        structure_volume = av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure;
        
        % Apply hemisphere filtering if specified
        if ~isempty(regionPlotLoc) && regionPlotLoc(iRegion) ~= 0
            % Allen CCF is in (AP, DV, ML) order when loaded
            % After permute for isosurface it becomes (ML, AP, DV)
            [ap_size, dv_size, ml_size] = size(structure_volume);
            ml_midline = round(ml_size / 2);
            
            if regionPlotLoc(iRegion) == -1
                % Keep only left hemisphere (ml < midline)
                structure_volume(:, :, ml_midline:end) = 0;
            elseif regionPlotLoc(iRegion) == 1
                % Keep only right hemisphere (ml > midline)
                structure_volume(:, :, 1:ml_midline) = 0;
            end
        end
        
        % Create isosurface from filtered volume
        structure_3d = isosurface(permute(structure_volume, [3, 1, 2]), 0);
        
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
    % Generate colors like MATLAB's lines() but with more variety
    numColors = max(23, length(theseAnimals));
    
    % Start with MATLAB's default lines colors and extend
    if numColors <= 7
        mouseColors = lines(numColors);
    else
        % Get the base lines colors
        baseColors = lines(7);
        
        % Generate additional colors with maximum separation like lines() does
        additionalColors = [];
        numAdditional = numColors - 7;
        
        % Create hues with maximum separation (like lines() algorithm)
        % Start with offsets to avoid overlapping with base colors
        hueOffsets = [0.15, 0.45, 0.75, 0.05, 0.35, 0.65, 0.95, 0.25, 0.55, 0.85, ...
                     0.10, 0.40, 0.70, 0.20, 0.50, 0.80, 0.30, 0.60, 0.90, 0.00];
        
        for i = 1:numAdditional
            % Use predefined hue offsets for maximum separation
            hue = hueOffsets(mod(i-1, length(hueOffsets)) + 1);
            
            % Vary saturation and brightness to maintain distinction
            sat = 0.7 + 0.2 * sin(i * pi/4); % Vary saturation between 0.7-0.9
            val = 0.75 + 0.2 * cos(i * pi/5); % Vary brightness between 0.75-0.95
            
            % Convert HSV to RGB
            hsv_color = [hue, sat, val];
            rgb_color = hsv2rgb(hsv_color);
            
            % Special handling for mouse 15 (index 8 in additionalColors) - make it less fluorescent
            if i == 8  % This would be mouse 15 (7 base colors + 8th additional = mouse 15)
                % Replace with a more muted color (darker blue-gray)
                rgb_color = [0.4, 0.5, 0.7];
            end
            
            additionalColors = [additionalColors; rgb_color];
        end
        
        % Combine base colors with additional colors
        mouseColors = [baseColors; additionalColors];
        
        % Use only what we need
        mouseColors = mouseColors(1:numColors, :);
        mouseColors = mouseColors(1:length(theseAnimals), :);
    end
    
    for iAnimal = 1:size(theseAnimals, 2)
            %iAnimal = iAnimal + 1;
            % Load probe data
            outputDir = ['/Users/jf5479/Dropbox/Histology/' theseAnimals{iAnimal} '/brainReg/'];
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
                                % Use Allen CCF structure tree since trajectory_areas uses Allen IDs
                                struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
                                if ~isempty(struct_idx)
                                    struct_curr = st.id(struct_idx(1));
                                    if any(probe_ccf(iProbe).trajectory_areas == struct_curr)
                                        plotThisProbe = true;
                                        break;
                                    end
                                end
                                end
                            end
                        end
                    end
                    
                    if plotThisProbe
                        % Determine which points to use based on onlyROIProbes setting
                        if onlyROIProbes && isfield(probe_ccf(iProbe), 'trajectory_coords') && ~isempty(probe_ccf(iProbe).trajectory_coords) && ...
                           isfield(probe_ccf(iProbe), 'trajectory_areas') && ~isempty(probe_ccf(iProbe).trajectory_areas)
                            
                            % Use trajectory_coords for ROI filtering (these match trajectory_areas)
                            % Get ROI structure IDs
                            roiStructureIds = [];
                            for iRegion = 1:length(regionsNames)
                                struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
                                if ~isempty(struct_idx)
                                    roiStructureIds = [roiStructureIds, st.id(struct_idx(1))];
                                end
                            end
                            
                            % Filter trajectory_coords to only those in ROI structures
                            if ~isempty(roiStructureIds)
                                trajectory_coords = probe_ccf(iProbe).trajectory_coords;
                                trajectory_areas = probe_ccf(iProbe).trajectory_areas;
                                
                                % Find which trajectory points are in ROI
                                pointsInROI = false(size(trajectory_coords, 1), 1);
                                for iPoint = 1:length(trajectory_areas)
                                    if any(roiStructureIds == trajectory_areas(iPoint))
                                        pointsInROI(iPoint) = true;
                                    end
                                end
                                
                                % Use filtered trajectory_coords if any ROI points found
                                if any(pointsInROI)
                                    % trajectory_coords are in (AP, DV, ML) order
                                    % but we need (AP, ML, DV) for plotting
                                    filtered_coords = trajectory_coords(pointsInROI, :);
                                    % Swap columns 2 and 3 to convert from (AP, DV, ML) to (AP, ML, DV)
                                    thesePoints = filtered_coords(:, [1, 3, 2]) * 2.5;
                                else
                                    continue; % Skip if no ROI points
                                end
                            else
                                continue; % Skip if no ROI structures defined
                            end
                        else
                            % Use the original sparse points (standard behavior)
                            thesePoints = probe_ccf(curr_probe).points * 2.5;
                        end
                        
                        % Check if we still have points after ROI filtering
                        if isempty(thesePoints)
                            continue; % Skip if no points remain after filtering
                        end
                        
                        % Detect multi-shank probe by looking for gaps in ML coordinates
                        % probe_ccf.points need to be rearranged for plotting
                        % Plot3 expects (AP, ML, DV) based on the axis limits
                        % Sort points by ML coordinate (2nd column) to find gaps
                        [sorted_ml, sort_idx] = sort(thesePoints(:, 2));
                        ml_diff = diff(sorted_ml);
                        
                        % Find large gaps (> 200um scaled) indicating separate shanks
                        shank_gap_threshold = 200; % 200um scaled gap between shanks
                        shank_boundaries = find(ml_diff > shank_gap_threshold);
                        
                        % Split points into shanks
                        if ~isempty(shank_boundaries)
                            % Multi-shank probe detected
                            n_shanks = length(shank_boundaries) + 1;
                            shank_points = cell(n_shanks, 1);
                            
                            % Assign points to shanks
                            shank_starts = [1; shank_boundaries + 1];
                            shank_ends = [shank_boundaries; length(sorted_ml)];
                            
                            for iShank = 1:n_shanks
                                shank_idx = sort_idx(shank_starts(iShank):shank_ends(iShank));
                                shank_points{iShank} = thesePoints(shank_idx, :);
                            end
                            
                        else
                            % Single shank probe
                            n_shanks = 1;
                            shank_points = {thesePoints};
                        end
                        
                        % Process each shank separately
                        for iShank = 1:n_shanks
                            theseShankPoints = shank_points{iShank};
                            
                            % Collect ALL hemisphere preferences for regions this probe passes through
                            targetHemispheres = [];
                            
                            if ~isempty(regionPlotLoc) && any(regionPlotLoc ~= 0)
                                % Check which regions this probe passes through
                                probePassesThroughROI = false;
                                if isfield(probe_ccf(iProbe), 'trajectory_areas') && ~isempty(probe_ccf(iProbe).trajectory_areas)
                                    for iRegion = 1:length(regionsNames)
                                        % Use Allen CCF structure tree since trajectory_areas uses Allen IDs
                                        struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
                                        if ~isempty(struct_idx)
                                            struct_curr = st.id(struct_idx(1));
                                            if any(probe_ccf(iProbe).trajectory_areas == struct_curr)
                                                probePassesThroughROI = true;
                                                if regionPlotLoc(iRegion) ~= 0
                                                    % Collect this region's hemisphere preference
                                                    targetHemispheres = [targetHemispheres, regionPlotLoc(iRegion)];
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                % Get unique hemisphere preferences
                                targetHemispheres = unique(targetHemispheres);
                                
                                % If onlyROIProbes is true and probe doesn't pass through ROI, skip it
                                if onlyROIProbes && ~probePassesThroughROI
                                    continue; % Skip this shank
                                end
                                
                                % If no hemisphere preferences found, use default
                                if isempty(targetHemispheres)
                                    targetHemispheres = 0; % Plot in actual position
                                end
                            else
                                targetHemispheres = 0; % Plot in actual position
                            end
                        
                            % Plot shank for each target hemisphere
                            for iHemisphere = 1:length(targetHemispheres)
                                currentTargetHemisphere = targetHemispheres(iHemisphere);
                                
                                % Use the shank points (already filtered if onlyROIProbes=true)
                                plotPoints = theseShankPoints;
                                
                                % Mirror points if needed
                                if currentTargetHemisphere ~= 0
                                % probe_ccf points are in CCF coordinates (AP, DV, ML)
                                % After scaling by 2.5, plotPoints are still (AP, DV, ML)
                                % In plot3(x, y, z), this maps to plot3(AP, DV, ML)
                                [ap_max, dv_max, ml_max] = size(tv);
                                
                                % Use the same midline as Witten lab code (570 after scaling)
                                ml_midline_scaled = 570; % This is the empirical midline value
                                
                                % Debug output - ML is 2nd coordinate
                                mean_ml = mean(plotPoints(:, 2));
                                is_left = mean_ml < ml_midline_scaled;
                                if is_left
                                    current_side = 'left';
                                else
                                    current_side = 'right';
                                end
                                
                                % Mirror ML coordinate (2nd coordinate) using correct formula
                                if currentTargetHemisphere == -1
                                    % Want probe on left hemisphere
                                    if ~is_left
                                        % Probe is on right, mirror to left
                                        plotPoints(:, 2) = 2*ml_midline_scaled - plotPoints(:, 2);
                                    end
                                elseif currentTargetHemisphere == 1
                                    % Want probe on right hemisphere  
                                    if is_left
                                        % Probe is on left, mirror to right
                                        plotPoints(:, 2) = 2*ml_midline_scaled - plotPoints(:, 2);
                                    end
                                end
                                
                            end
                            
                            % Plot probe points if requested
                            if showPoints
                                plot3(plotPoints(:, 1), ...
                                    plotPoints(:, 2), ...
                                    plotPoints(:, 3), ...
                                    '.', 'color', animalColor, 'MarkerSize', 20);
                            end
                            
                            % Fit curve through points
                            if useBezierFit && size(plotPoints, 1) >= 3
                                % Use Bezier curve fitting
                                % Sort points by one dimension (e.g., DV) to get proper curve order
                                [~, sort_idx] = sort(plotPoints(:, 3)); % Sort by DV (z-axis)
                                sorted_points = plotPoints(sort_idx, :);
                                
                                % Create Bezier curve
                                t = linspace(0, 1, 1000);
                                bezier_curve_points = bezier_curve(t, sorted_points);
                                
                                % Mirror bezier curve if needed (in case it extends beyond probe points)
                                if currentTargetHemisphere ~= 0
                                    ml_midline_scaled_bezier = 570; % Same empirical midline
                                    mean_ml_bezier = mean(bezier_curve_points(:, 2)); % ML is 2nd coordinate
                                    is_left_bezier = mean_ml_bezier < ml_midline_scaled_bezier;
                                    
                                    if currentTargetHemisphere == -1 && ~is_left_bezier
                                        % Want on left but is on right
                                        bezier_curve_points(:, 2) = 2*ml_midline_scaled_bezier - bezier_curve_points(:, 2);
                                    elseif currentTargetHemisphere == 1 && is_left_bezier
                                        % Want on right but is on left
                                        bezier_curve_points(:, 2) = 2*ml_midline_scaled_bezier - bezier_curve_points(:, 2);
                                    end
                                end
                                
                                % Plot Bezier curve
                                plot3(bezier_curve_points(:, 1), bezier_curve_points(:, 2), bezier_curve_points(:, 3), ...
                                    'color', animalColor, 'linewidth', 2);
                            else
                                % Use linear fit but only within the ROI segment bounds
                                if size(plotPoints, 1) >= 2
                                    r0 = mean(plotPoints, 1);
                                    xyz = bsxfun(@minus, plotPoints, r0);
                                    [~, ~, V] = svd(xyz, 0);
                                    histology_probe_direction = V(:, 1);
                                    
                                    % Make sure the direction goes down in DV - flip if it's going up
                                    if histology_probe_direction(3) < 0
                                        histology_probe_direction = -histology_probe_direction;
                                    end

                                    % Limit line to the extent of the ROI points instead of extending infinitely
                                    minPoint = min(plotPoints, [], 1);
                                    maxPoint = max(plotPoints, [], 1);
                                    
                                    % Project endpoints onto the probe direction
                                    range = max(sqrt(sum((plotPoints - r0).^2, 2))) * 1.1; % Extend slightly beyond data
                                    line_eval = [-range, range];
                                    probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0);
                                    
                                    % Plot linear fit
                                    plot3(probe_fit_line(:, 1), probe_fit_line(:, 2), probe_fit_line(:, 3), ...
                                        'color', animalColor, 'linewidth', 2);
                                else
                                    % Just plot the points if less than 2 points
                                    plot3(plotPoints(:, 1), plotPoints(:, 2), plotPoints(:, 3), ...
                                        'color', animalColor, 'linewidth', 2);
                                end
                            end
                            end % end hemisphere loop
                        end % end shank loop
                    end % end if plotThisProbe
                    
                    % Count probes per region (regardless of plotting)
                    if isfield(probe_ccf(iProbe), 'trajectory_areas')
                        for iRegion = 1:length(regionsNames)
                            % Use Allen CCF structure tree since trajectory_areas uses Allen IDs
                            struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
                            if ~isempty(struct_idx)
                                struct_curr = st.id(struct_idx(1));
                                if any(probe_ccf(iProbe).trajectory_areas == struct_curr)
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
            end
    end % end animal loop
 % end animalsType loop

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
t = title(sprintf('Probe tracks for %d mice in %s', length(theseAnimals), strjoin(regionsNames, ', ')));
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
    legendEntries{end+1} = sprintf('Mouse %d', iMouse);
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
    plotProbeRegions(theseAnimals, mouseColors, onlyROIProbes, regionsNames, st_br, allenAtlasPath, st);
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

function roiPoints = filterPointsToROI(points, regionsNames, st_br, probe_ccf_struct, st)
% Filter probe points to only those within the regions of interest
% Uses the existing trajectory_areas data instead of looking up coordinates
% points: Nx3 array of probe points (already scaled by 2.5)
% probe_ccf_struct: the probe_ccf structure containing trajectory_areas
% Returns only the points that fall within any of the ROI structures

roiPoints = [];

if isempty(points) || ~isfield(probe_ccf_struct, 'trajectory_areas') || isempty(probe_ccf_struct.trajectory_areas)
    return;
end

% Get ROI structure IDs from Allen CCF structure tree (st)
% trajectory_areas uses Allen CCF IDs, not brainglobe IDs
roiStructureIds = [];
for iRegion = 1:length(regionsNames)
    struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
    if ~isempty(struct_idx)
        struct_id = st.id(struct_idx(1));
        roiStructureIds = [roiStructureIds, struct_id];
    end
end

if isempty(roiStructureIds)
    return;
end

% Check each point using trajectory_areas data
trajectory_areas = probe_ccf_struct.trajectory_areas;

% We need to map from the input points to the corresponding trajectory_areas
% The challenge is that points might be a subset (e.g., shank points) of the full trajectory
% For simplicity, we'll assume points correspond to the first N trajectory areas
% This works for most cases where points are extracted sequentially

pointsInROI = false(size(points, 1), 1);
numPointsToCheck = min(size(points, 1), length(trajectory_areas));

for iPoint = 1:numPointsToCheck
    structureId = trajectory_areas(iPoint);
    if any(roiStructureIds == structureId)
        pointsInROI(iPoint) = true;
    end
end

% Return only points that are in ROI structures
if any(pointsInROI)
    roiPoints = points(pointsInROI, :);
else
    roiPoints = [];
end

end

function plotProbeRegions(theseAnimals, mouseColors, onlyROIProbes, regionsNames, st_br, allenAtlasPath, st)
% Plot brain regions that each probe passes through
% Load Allen CCF colormap
cmap_filename = [allenAtlasPath, filesep, 'allen_ccf_colormap_2017.mat'];
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
    % Load probe data using the same path as main function
    outputDir = ['/Users/jf5479/Dropbox/Histology/' theseAnimals{iAnimal} '/brainReg/'];
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
                    % Use Allen CCF structure tree since trajectory_areas uses Allen IDs
                    struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
                    if ~isempty(struct_idx)
                        struct_curr = st.id(struct_idx(1));
                        if any(probe_ccf(iProbe).trajectory_areas == struct_curr)
                            includeThisProbe = true;
                            break;
                        end
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
            % Use Allen CCF structure tree since trajectory_areas uses Allen IDs
            matching_idx = find(st.id == area_id);
            if ~isempty(matching_idx)
                trajectory_area_labels{iArea} = st.acronym{matching_idx(1)};
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
%sgtitle('Brain Regions Traversed by Each Probe', 'FontSize', 14, 'Color', 'w');

end

