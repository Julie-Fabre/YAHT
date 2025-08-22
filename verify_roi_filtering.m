% Verify that ROI filtering is actually working
clear all;
cd('/Users/jf5479/Dropbox/MATLAB/onPaths/YAHT');
cl_myPaths;

fprintf('=== Verifying ROI Filtering Implementation ===\n');

% Load structure trees
st = ya_loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']);
regionsNames = {'CP', 'GPe', 'SNr'};

% Get ROI structure IDs
roiStructureIds = [];
for iRegion = 1:length(regionsNames)
    struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
    if ~isempty(struct_idx)
        roiStructureIds = [roiStructureIds, st.id(struct_idx(1))];
        fprintf('%s: ID = %d\n', regionsNames{iRegion}, st.id(struct_idx(1)));
    end
end

% Test with one animal
animal = 'JF107';
outputDir = ['/Users/jf5479/Dropbox/Histology/' animal '/brainReg/'];
probe_ccf_location = [outputDir, 'probe_ccf.mat'];
load(probe_ccf_location);

fprintf('\n=== Analyzing %s ===\n', animal);

for iProbe = 1:min(3, length(probe_ccf)) % Check first 3 probes
    if isfield(probe_ccf(iProbe), 'points') && ~isempty(probe_ccf(iProbe).points) && ...
       isfield(probe_ccf(iProbe), 'trajectory_areas') && ~isempty(probe_ccf(iProbe).trajectory_areas)
        
        fprintf('\nProbe %d:\n', iProbe);
        
        % Original points
        origPoints = probe_ccf(iProbe).points;
        fprintf('  Original: %d points\n', size(origPoints, 1));
        
        % Count points in ROI
        trajectory_areas = probe_ccf(iProbe).trajectory_areas;
        numPoints = min(size(origPoints, 1), length(trajectory_areas));
        pointsInROI = false(numPoints, 1);
        
        for iPoint = 1:numPoints
            if any(roiStructureIds == trajectory_areas(iPoint))
                pointsInROI(iPoint) = true;
            end
        end
        
        numROIPoints = sum(pointsInROI);
        fprintf('  In ROI: %d points (%.1f%%)\n', numROIPoints, 100*numROIPoints/numPoints);
        
        % Show which regions
        for iRegion = 1:length(regionsNames)
            struct_idx = find(strcmp(st.acronym, regionsNames{iRegion}));
            if ~isempty(struct_idx)
                struct_id = st.id(struct_idx(1));
                if any(trajectory_areas == struct_id)
                    fprintf('    Passes through %s\n', regionsNames{iRegion});
                end
            end
        end
    end
end

fprintf('\n=== Visual Test ===\n');
fprintf('Calling function with onlyROIProbes=true...\n');
fprintf('You should see ONLY probe segments within CP, GPe, SNr\n');
fprintf('NOT the full dorsal-ventral extent of probes\n\n');

ya_plotAllProbeTracksInROIs_JF({'JF107'}, regionsNames, [], true, true, false, false, [], [], [], [-1, 1, -1]);

fprintf('\n✓ If probes appear truncated to only ROI regions, filtering is working!\n');
fprintf('✗ If probes show full dorsal-ventral extent, filtering is NOT working.\n');