% Minimal test of function initialization
clear all;

fprintf('=== Testing function initialization ===\n');

% Test parameters
theseAnimals = {'JF107'};
regionsOfInterest = {'CP', 'GPe', 'SNr'};
plotHemisphere = [-1, 1, -1];

try
    fprintf('Step 1: Setting up parameters...\n');
    
    % Replicate the parameter setup from the function
    onlyROIProbes = 1;
    showPoints = 1;
    useBezierFit = 0;
    showRegionPlot = 0;  % Disable for now
    regionColors = 'allen';
    blackBackground = 0;
    thickBrainLines = 0.5;
    
    % Initialize
    animalsType = {'Naive'};
    regionsNames = regionsOfInterest;
    regions = regionsOfInterest;
    
    % Set up hemisphere plotting
    if ~isempty(plotHemisphere)
        if length(plotHemisphere) ~= length(regionsOfInterest)
            error('plotHemisphere must have the same length as regionsOfInterest');
        end
        regionPlotLoc = plotHemisphere;
    else
        regionPlotLoc = zeros(1, length(regionsNames));
    end
    
    fprintf('✓ Parameters set up successfully\n');
    
    fprintf('Step 2: Loading paths...\n');
    cl_myPaths;
    fprintf('✓ Paths loaded\n');
    
    fprintf('Step 3: Loading atlas data...\n');
    tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']);
    fprintf('✓ Template volume loaded: %dx%dx%d\n', size(tv,1), size(tv,2), size(tv,3));
    
    av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']);
    fprintf('✓ Annotation volume loaded\n');
    
    st = ya_loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']);
    fprintf('✓ Allen structure tree loaded: %d structures\n', length(st.id));
    
    % Load brainglobe atlas
    atlasType = 'allen';
    atlasSpecies = 'mouse';
    atlasResolution_um = 25;
    atlasLocation = dir([brainglobeLocation, atlasType, '_', atlasSpecies, '_', num2str(atlasResolution_um), 'um*']);
    [~, ~, st_br, ~] = ya_loadAllenAtlas([atlasLocation.folder, filesep, atlasLocation.name]);
    fprintf('✓ Brainglobe structure tree loaded: %d structures\n', length(st_br.id));
    
    fprintf('Step 4: Testing structure ID lookup...\n');
    for iRegion = 1:length(regionsNames)
        % Allen CCF
        struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
        if ~isempty(struct_idx)
            allen_id = st.id(struct_idx(1));
            fprintf('  %s: Allen ID = %d\n', regionsNames{iRegion}, allen_id);
        else
            fprintf('  %s: NOT FOUND in Allen tree!\n', regionsNames{iRegion});
        end
        
        % Brainglobe
        struct_idx_br = find(strcmp(st_br.acronym, regionsNames{iRegion}));
        if ~isempty(struct_idx_br)
            bg_id = st_br.id(struct_idx_br(1));
            fprintf('  %s: Brainglobe ID = %d\n', regionsNames{iRegion}, bg_id);
        else
            fprintf('  %s: NOT FOUND in Brainglobe tree!\n', regionsNames{iRegion});
        end
    end
    
    fprintf('Step 5: Testing probe file loading...\n');
    animal = theseAnimals{1};
    outputDir = ['/Users/jf5479/Dropbox/Histology/' animal '/brainReg/'];
    probe_ccf_location = [outputDir, 'probe_ccf.mat'];
    
    if exist(probe_ccf_location, 'file')
        load(probe_ccf_location);
        fprintf('✓ Loaded %s with %d probes\n', animal, length(probe_ccf));
        
        % Test ROI detection on first probe with data
        for iProbe = 1:length(probe_ccf)
            if isfield(probe_ccf(iProbe), 'trajectory_areas') && ~isempty(probe_ccf(iProbe).trajectory_areas)
                fprintf('  Probe %d: %d trajectory areas\n', iProbe, length(probe_ccf(iProbe).trajectory_areas));
                
                % Check against Allen IDs
                for iRegion = 1:length(regionsNames)
                    struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
                    if ~isempty(struct_idx)
                        struct_curr = st.id(struct_idx(1));
                        if any(probe_ccf(iProbe).trajectory_areas == struct_curr)
                            fprintf('    ✓ Passes through %s (Allen ID %d)\n', regionsNames{iRegion}, struct_curr);
                        end
                    end
                end
                break; % Just test first probe with data
            end
        end
    else
        fprintf('✗ Probe file not found: %s\n', probe_ccf_location);
    end
    
    fprintf('\n✓ All initialization steps completed successfully!\n');
    
catch e
    fprintf('\n✗ ERROR at step: %s\n', e.message);
    if ~isempty(e.stack)
        fprintf('Line %d in %s\n', e.stack(1).line, e.stack(1).name);
    end
end