% Debug script to understand probe_ccf coordinate system
clear all;
cl_myPaths;

% Load atlas dimensions
tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']);
[ap_max, dv_max, ml_max] = size(tv);
fprintf('Atlas dimensions (AP, DV, ML): [%d, %d, %d]\n', ap_max, dv_max, ml_max);

% Load a probe file
probe_file = '/Users/jf5479/Dropbox/Histology/JF109/brainReg/probe_ccf.mat';
load(probe_file);

% Check first probe with points
for i = 1:length(probe_ccf)
    if isfield(probe_ccf(i), 'points') && ~isempty(probe_ccf(i).points)
        points = probe_ccf(i).points;
        fprintf('\nProbe %d raw coordinates:\n', i);
        fprintf('  Col 1 range: [%.1f, %.1f] (compare to AP max=%d)\n', min(points(:,1)), max(points(:,1)), ap_max);
        fprintf('  Col 2 range: [%.1f, %.1f] (compare to DV max=%d)\n', min(points(:,2)), max(points(:,2)), dv_max);
        fprintf('  Col 3 range: [%.1f, %.1f] (compare to ML max=%d)\n', min(points(:,3)), max(points(:,3)), ml_max);
        
        % The probe should span a large range in DV (dorsal-ventral)
        % And have limited range in AP and ML
        col1_span = max(points(:,1)) - min(points(:,1));
        col2_span = max(points(:,2)) - min(points(:,2));
        col3_span = max(points(:,3)) - min(points(:,3));
        
        fprintf('\n  Spans:\n');
        fprintf('    Col 1 span: %.1f\n', col1_span);
        fprintf('    Col 2 span: %.1f\n', col2_span);
        fprintf('    Col 3 span: %.1f\n', col3_span);
        
        % Largest span should be DV (probes go dorsal to ventral)
        [~, dv_col] = max([col1_span, col2_span, col3_span]);
        fprintf('\n  Column %d has largest span - likely DV axis\n', dv_col);
        
        % Check hemisphere
        if col3_span < 100 % Narrow in ML
            ml_mean = mean(points(:,3));
            if ml_mean < ml_max/2
                fprintf('  Based on col 3: LEFT hemisphere (mean=%.1f < midline=%.1f)\n', ml_mean, ml_max/2);
            else
                fprintf('  Based on col 3: RIGHT hemisphere (mean=%.1f > midline=%.1f)\n', ml_mean, ml_max/2);
            end
        end
        
        if col2_span < 100 % Narrow in what might be ML
            ml_mean = mean(points(:,2));
            if ml_mean < dv_max/2
                fprintf('  If col 2 is ML: LEFT hemisphere (mean=%.1f < midline=%.1f)\n', ml_mean, dv_max/2);
            else
                fprintf('  If col 2 is ML: RIGHT hemisphere (mean=%.1f > midline=%.1f)\n', ml_mean, dv_max/2);
            end
        end
        
        break; % Just check first probe with points
    end
end

fprintf('\n=== CONCLUSION ===\n');
fprintf('probe_ccf.points columns are most likely:\n');
fprintf('  Column 1: AP (anterior-posterior)\n');
fprintf('  Column 2: DV (dorsal-ventral) - largest span for vertical probes\n');
fprintf('  Column 3: ML (medial-lateral)\n');
fprintf('\nBut the plot3 visualization might use a different mapping!\n');