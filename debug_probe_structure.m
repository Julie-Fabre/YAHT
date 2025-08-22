% Debug probe structure to understand points vs trajectory_areas
clear all;
cd('/Users/jf5479/Dropbox/MATLAB/onPaths/YAHT');
cl_myPaths;

% Load one probe file
animal = 'JF107';
outputDir = ['/Users/jf5479/Dropbox/Histology/' animal '/brainReg/'];
probe_ccf_location = [outputDir, 'probe_ccf.mat'];
load(probe_ccf_location);

fprintf('=== Debugging Probe Structure ===\n');

% Check probe 3 which passes through CP
iProbe = 3;
fprintf('\nProbe %d analysis:\n', iProbe);

if isfield(probe_ccf(iProbe), 'points')
    fprintf('  points: %d x %d\n', size(probe_ccf(iProbe).points));
end

if isfield(probe_ccf(iProbe), 'trajectory_areas')
    fprintf('  trajectory_areas: %d elements\n', length(probe_ccf(iProbe).trajectory_areas));
end

if isfield(probe_ccf(iProbe), 'trajectory_coords')
    fprintf('  trajectory_coords: %d x %d\n', size(probe_ccf(iProbe).trajectory_coords));
end

% The key insight: points might not correspond 1:1 with trajectory_areas!
% trajectory_areas likely corresponds to trajectory_coords, not points

fprintf('\n=== Key Finding ===\n');
fprintf('points (%d) and trajectory_areas (%d) have different lengths!\n', ...
    size(probe_ccf(iProbe).points, 1), length(probe_ccf(iProbe).trajectory_areas));
fprintf('trajectory_coords (%d) matches trajectory_areas (%d)\n', ...
    size(probe_ccf(iProbe).trajectory_coords, 1), length(probe_ccf(iProbe).trajectory_areas));

% Check what trajectory_coords represents
fprintf('\n=== Understanding the data ===\n');
fprintf('probe_ccf.points: Likely sparse probe track points from manual annotation\n');
fprintf('probe_ccf.trajectory_coords: Dense trajectory points (one per trajectory_area)\n');
fprintf('probe_ccf.trajectory_areas: Brain region ID for each trajectory_coord\n');

% Verify CP is in trajectory_areas
st = ya_loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']);
cp_id = st.id(find(strcmp(st.acronym, 'CP'), 1));
fprintf('\nCP (ID=%d) appears %d times in trajectory_areas\n', ...
    cp_id, sum(probe_ccf(iProbe).trajectory_areas == cp_id));