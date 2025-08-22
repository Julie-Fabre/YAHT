% Comprehensive debug to understand the probe detection issue
clear all;
cl_myPaths;

% Use exact same parameters as user's call
theseAnimals = {'JF107', 'JF093', 'JF091', 'JF104', 'JF105', 'JF118', 'JF119', 'JF120', 'JF121'};
regionsNames = {'CP', 'GPe', 'SNr'};
onlyROIProbes = true; % Assuming this is true based on user's issue
plotHemisphere = [-1, 1, -1];

% Load structure trees
fprintf('Loading structure trees...\n');
st = ya_loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']);
atlasType = 'allen';
atlasSpecies = 'mouse';
atlasResolution_um = 25;
atlasLocation = dir([brainglobeLocation, atlasType, '_', atlasSpecies, '_', num2str(atlasResolution_um), 'um*']);
[~, ~, st_br, ~] = ya_loadAllenAtlas([atlasLocation.folder, filesep, atlasLocation.name]);

% Show structure IDs for our regions
fprintf('\n=== Structure IDs ===\n');
for iRegion = 1:length(regionsNames)
    % Allen CCF
    struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
    allen_id = st.id(struct_idx(1));
    
    % Brainglobe  
    struct_idx_br = find(strcmp(st_br.acronym, regionsNames{iRegion}));
    bg_id = st_br.id(struct_idx_br(1));
    
    fprintf('%s: Allen=%d, Brainglobe=%d\n', regionsNames{iRegion}, allen_id, bg_id);
end

% Test each animal
totalProbesFound = 0;
totalProbesWithTrajectory = 0;
totalProbesPassingROI = 0;

for iAnimal = 1:length(theseAnimals)
    animal = theseAnimals{iAnimal};
    fprintf('\n=== Animal %s ===\n', animal);
    
    % Try to load probe file
    probe_ccf_location = cl_cortexlab_filename(animal, [], [], 'histo');
    fprintf('Probe file: %s\n', probe_ccf_location);
    
    if ~exist(probe_ccf_location, 'file')
        fprintf('  FILE NOT FOUND!\n');
        continue;
    end
    
    load(probe_ccf_location);
    fprintf('  Loaded %d probes\n', length(probe_ccf));
    
    for iProbe = 1:length(probe_ccf)
        totalProbesFound = totalProbesFound + 1;
        
        % Check basic data
        hasPoints = isfield(probe_ccf(iProbe), 'points') && ~isempty(probe_ccf(iProbe).points);
        hasTrajectory = isfield(probe_ccf(iProbe), 'trajectory_areas') && ~isempty(probe_ccf(iProbe).trajectory_areas);
        
        if hasPoints && hasTrajectory
            totalProbesWithTrajectory = totalProbesWithTrajectory + 1;
            
            % Check trajectory areas
            trajectory_areas = probe_ccf(iProbe).trajectory_areas;
            unique_areas = unique(trajectory_areas);
            
            fprintf('    Probe %d: %d points, %d trajectory areas\n', iProbe, size(probe_ccf(iProbe).points, 1), length(trajectory_areas));
            fprintf('      Sample areas: %s\n', sprintf('%d ', unique_areas(1:min(5, length(unique_areas)))));
            
            % Check against Allen CCF structure IDs (correct method)
            passesROI_allen = false;
            for iRegion = 1:length(regionsNames)
                struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
                if ~isempty(struct_idx)
                    struct_id = st.id(struct_idx(1));
                    if any(trajectory_areas == struct_id)
                        fprintf('      PASSES %s (Allen ID %d)\n', regionsNames{iRegion}, struct_id);
                        passesROI_allen = true;
                    end
                end
            end
            
            % Check against Brainglobe structure IDs (wrong method - for comparison)
            passesROI_bg = false;
            for iRegion = 1:length(regionsNames)
                struct_idx_br = find(strcmp(st_br.acronym, regionsNames{iRegion}));
                if ~isempty(struct_idx_br)
                    struct_id = st_br.id(struct_idx_br(1));
                    if any(trajectory_areas == struct_id)
                        fprintf('      PASSES %s (Brainglobe ID %d) - OLD METHOD\n', regionsNames{iRegion}, struct_id);
                        passesROI_bg = true;
                    end
                end
            end
            
            if passesROI_allen
                totalProbesPassingROI = totalProbesPassingROI + 1;
                fprintf('      --> This probe SHOULD be plotted\n');
            else
                fprintf('      --> This probe should be skipped\n');
            end
            
        else
            fprintf('    Probe %d: missing data (points=%d, trajectory=%d)\n', iProbe, hasPoints, hasTrajectory);
        end
    end
end

fprintf('\n=== SUMMARY ===\n');
fprintf('Total probes found: %d\n', totalProbesFound);
fprintf('Probes with valid data: %d\n', totalProbesWithTrajectory);
fprintf('Probes passing ROI (Allen method): %d\n', totalProbesPassingROI);
fprintf('Expected result: %d probes should be plotted\n', totalProbesPassingROI);