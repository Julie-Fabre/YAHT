% Debug coordinate systems for points vs trajectory_coords
clear all;
cd('/Users/jf5479/Dropbox/MATLAB/onPaths/YAHT');
cl_myPaths;

% Load probe data
animal = 'JF107';
outputDir = ['/Users/jf5479/Dropbox/Histology/' animal '/brainReg/'];
probe_ccf_location = [outputDir, 'probe_ccf.mat'];
load(probe_ccf_location);

% Load atlas to get dimensions
tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']);
[ap_max, dv_max, ml_max] = size(tv);
fprintf('Atlas dimensions (AP, DV, ML): [%d, %d, %d]\n', ap_max, dv_max, ml_max);

% Find a probe with both points and trajectory_coords
for iProbe = 1:length(probe_ccf)
    if isfield(probe_ccf(iProbe), 'points') && ~isempty(probe_ccf(iProbe).points) && ...
       isfield(probe_ccf(iProbe), 'trajectory_coords') && ~isempty(probe_ccf(iProbe).trajectory_coords)
        
        fprintf('\n=== Probe %d Coordinate Analysis ===\n', iProbe);
        
        points = probe_ccf(iProbe).points;
        trajectory_coords = probe_ccf(iProbe).trajectory_coords;
        
        fprintf('\nOriginal points (%d total):\n', size(points, 1));
        fprintf('  Column 1 range: [%.1f, %.1f] span=%.1f\n', min(points(:,1)), max(points(:,1)), range(points(:,1)));
        fprintf('  Column 2 range: [%.1f, %.1f] span=%.1f\n', min(points(:,2)), max(points(:,2)), range(points(:,2)));
        fprintf('  Column 3 range: [%.1f, %.1f] span=%.1f\n', min(points(:,3)), max(points(:,3)), range(points(:,3)));
        
        fprintf('\nTrajectory_coords (%d total):\n', size(trajectory_coords, 1));
        fprintf('  Column 1 range: [%.1f, %.1f] span=%.1f\n', min(trajectory_coords(:,1)), max(trajectory_coords(:,1)), range(trajectory_coords(:,1)));
        fprintf('  Column 2 range: [%.1f, %.1f] span=%.1f\n', min(trajectory_coords(:,2)), max(trajectory_coords(:,2)), range(trajectory_coords(:,2)));
        fprintf('  Column 3 range: [%.1f, %.1f] span=%.1f\n', min(trajectory_coords(:,3)), max(trajectory_coords(:,3)), range(trajectory_coords(:,3)));
        
        % Identify which column is which based on span
        fprintf('\n=== Coordinate Identification ===\n');
        
        % For a vertical probe, DV should have the largest span
        [~, points_dv_col] = max([range(points(:,1)), range(points(:,2)), range(points(:,3))]);
        [~, traj_dv_col] = max([range(trajectory_coords(:,1)), range(trajectory_coords(:,2)), range(trajectory_coords(:,3))]);
        
        fprintf('Points: Column %d has largest span (likely DV)\n', points_dv_col);
        fprintf('Trajectory: Column %d has largest span (likely DV)\n', traj_dv_col);
        
        % The plotting function uses plot3(AP, ML, DV)
        % Based on the Witten lab code, after scaling points are used as:
        % plot3(points(:,1), points(:,2), points(:,3))
        % Which suggests points are already in (AP, ML, DV) order
        
        fprintf('\n=== Expected coordinate mapping ===\n');
        fprintf('For points: Column 1=AP, Column 2=ML, Column 3=DV\n');
        fprintf('For trajectory_coords: Need to determine...\n');
        
        % Check if trajectory_coords matches points orientation
        if size(points, 1) > 1 && size(trajectory_coords, 1) > 1
            % Sample first and last points
            points_start = points(1, :);
            points_end = points(end, :);
            traj_start = trajectory_coords(1, :);
            traj_end = trajectory_coords(end, :);
            
            fprintf('\nFirst point comparison:\n');
            fprintf('  points:     [%.1f, %.1f, %.1f]\n', points_start);
            fprintf('  trajectory: [%.1f, %.1f, %.1f]\n', traj_start);
            
            fprintf('\nLast point comparison:\n');
            fprintf('  points:     [%.1f, %.1f, %.1f]\n', points_end);
            fprintf('  trajectory: [%.1f, %.1f, %.1f]\n', traj_end);
        end
        
        break; % Just analyze first valid probe
    end
end