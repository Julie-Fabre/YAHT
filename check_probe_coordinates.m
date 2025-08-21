% Check probe coordinate ranges to understand the coordinate system
clear all;

% Load Allen atlas to check dimensions
cl_myPaths;
tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']);
[ap_max, dv_max, ml_max] = size(tv);

fprintf('\n=== Allen CCF Template Volume Dimensions ===\n');
fprintf('AP dimension: %d voxels (%.1f mm)\n', ap_max, ap_max*0.01);
fprintf('DV dimension: %d voxels (%.1f mm)\n', dv_max, dv_max*0.01);
fprintf('ML dimension: %d voxels (%.1f mm)\n', ml_max, ml_max*0.01);
fprintf('Expected ML midline: %d voxels\n', ml_max/2);
fprintf('Expected ML midline scaled by 2.5: %.1f\n\n', ml_max/2*2.5);

% Check a few animals' probe coordinates
animals = {'JF109', 'JF110', 'JF059'};

for iAnimal = 1:length(animals)
    probe_file = ['/Users/jf5479/Dropbox/Histology/' animals{iAnimal} '/brainReg/probe_ccf.mat'];
    
    if exist(probe_file, 'file')
        fprintf('=== %s ===\n', animals{iAnimal});
        load(probe_file);
        
        for iProbe = 1:min(3, length(probe_ccf)) % Check first 3 probes
            if isfield(probe_ccf(iProbe), 'points') && ~isempty(probe_ccf(iProbe).points)
                points = probe_ccf(iProbe).points;
                
                % Check coordinate ranges (unscaled)
                ap_range = [min(points(:,1)), max(points(:,1))];
                dv_range = [min(points(:,2)), max(points(:,2))];
                ml_range = [min(points(:,3)), max(points(:,3))];
                
                fprintf('Probe %d:\n', iProbe);
                fprintf('  AP range: [%.1f, %.1f] voxels\n', ap_range);
                fprintf('  DV range: [%.1f, %.1f] voxels\n', dv_range);
                fprintf('  ML range: [%.1f, %.1f] voxels, mean=%.1f\n', ml_range, mean(points(:,3)));
                
                % Check if multi-shank by looking for ML gaps
                points_scaled = points * 2.5;
                [sorted_ml, ~] = sort(points_scaled(:,2));
                ml_diff = diff(sorted_ml);
                large_gaps = find(ml_diff > 200);
                if ~isempty(large_gaps)
                    fprintf('  Multi-shank detected! Gap size: %.1f\n', max(ml_diff));
                end
                
                % Determine hemisphere
                ml_mean = mean(points(:,3));
                if ml_mean < ml_max/2
                    hemisphere = 'LEFT';
                else
                    hemisphere = 'RIGHT';
                end
                fprintf('  Hemisphere: %s (ML mean %.1f vs midline %.1f)\n\n', hemisphere, ml_mean, ml_max/2);
            end
        end
    end
end

fprintf('\n=== COORDINATE SYSTEM ANALYSIS ===\n');
fprintf('If ML values are around 200-900, then:\n');
fprintf('  - ML axis goes from 0 (left) to ~1140 (right)\n');
fprintf('  - Midline is at ~570\n');
fprintf('  - Left hemisphere: ML < 570\n');
fprintf('  - Right hemisphere: ML > 570\n');