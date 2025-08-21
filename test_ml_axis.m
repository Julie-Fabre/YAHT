% Quick test to determine which axis is ML in the plotting
clear all;
cl_myPaths;

% Load template to get dimensions
tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']);
[ap_max, dv_max, ml_max] = size(tv);

fprintf('Template dimensions (AP, DV, ML): [%d, %d, %d]\n', ap_max, dv_max, ml_max);

% Load one probe to check
probe_file = '/Users/jf5479/Dropbox/Histology/JF109/brainReg/probe_ccf.mat';
load(probe_file);

if ~isempty(probe_ccf(1).points)
    points = probe_ccf(1).points;
    fprintf('\nProbe 1 coordinate ranges:\n');
    fprintf('Column 1 (AP?): [%.1f, %.1f]\n', min(points(:,1)), max(points(:,1)));
    fprintf('Column 2 (DV?): [%.1f, %.1f]\n', min(points(:,2)), max(points(:,2)));
    fprintf('Column 3 (ML?): [%.1f, %.1f]\n', min(points(:,3)), max(points(:,3)));
    
    % After scaling by 2.5
    points_scaled = points * 2.5;
    fprintf('\nAfter scaling by 2.5:\n');
    fprintf('Column 1: [%.1f, %.1f]\n', min(points_scaled(:,1)), max(points_scaled(:,1)));
    fprintf('Column 2: [%.1f, %.1f]\n', min(points_scaled(:,2)), max(points_scaled(:,2)));
    fprintf('Column 3: [%.1f, %.1f]\n', min(points_scaled(:,3)), max(points_scaled(:,3)));
    
    % The probe_ccf.points should be in CCF coordinates: (AP, DV, ML)
    % But in plot3, the order is plot3(x, y, z)
    % Looking at the brain visualization, we need to determine the mapping
    fprintf('\n=== ANALYSIS ===\n');
    fprintf('probe_ccf.points columns are likely: [AP, DV, ML]\n');
    fprintf('In plot3(points(:,1), points(:,2), points(:,3)):\n');
    fprintf('  x-axis (1st arg) = AP\n');
    fprintf('  y-axis (2nd arg) = DV\n');
    fprintf('  z-axis (3rd arg) = ML\n');
    fprintf('\nBut wait! Let me check if there was a permutation...\n');
end