function ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsOfInterest, patchBrain, onlyROIProbes)
% ya_plotAllProbeTracksInROIs_JF - Plot probe tracks for multiple animals with regions of interest
%
% Inputs:
%   theseAnimals - Cell array of animal names (e.g., {'JF058', 'JF059'})
%   regionsOfInterest - Cell array of region names to highlight (e.g., {'CP', 'GPe', 'SNr'})
%                       If empty or not provided, defaults to {'CP', 'GPe', 'GPi', 'STN', 'SNr'}
%   patchBrain - Boolean to use surface patch (1) or wire grid (0) for brain. Default: 0
%   onlyROIProbes - Boolean to plot only probes that pass through ROIs (1) or all probes (0). Default: 0
%
% Example:
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 1); % Only ROI probes
%   ya_plotAllProbeTracksInROIs_JF({'JF058', 'JF059'}, {'CP', 'SNr'}, 0, 0); % All probes

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
theseColors = lines(length(regionsNames));
theseColors = mat2cell(theseColors, ones(size(theseColors,1),1), 3);
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
        [~, brain_outline] = plotBrainGrid([], []);
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

            %thisAnimal = strcmp(recordingInfo.Mouse, qqqq{iAnimal});
            %ttP = (recordingInfo.HistologyProbe( thisAnimal & theseProbes));
            %thisProbe = recordingInfo.HistologyProbe(find(thisAnimal & theseProbes));
            % Get colors for this animal
            animalColors = lines(size(probe_ccf, 1));
            
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
                        
                        % Fit line through points
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
                        
                        % Plot probe points and fitted line
                        plot3(thesePoints(:, 1), ...
                            thesePoints(:, 2), ...
                            thesePoints(:, 3), ...
                            '.', 'color', animalColors(iProbe, :), 'MarkerSize', 20);
                        line(probe_fit_line(:, 1), probe_fit_line(:, 2), probe_fit_line(:, 3), ...
                            'color', animalColors(iProbe, :), 'linewidth', 2)
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
set(gcf, 'Color', 'w');

% Enable 3D rotation
h = rotate3d(gca);
h.Enable = 'on';

% Add title
title(sprintf('Probe tracks for %d animals in %s', length(theseAnimals), strjoin(regionsNames, ', ')));

end

