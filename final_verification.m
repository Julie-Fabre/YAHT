% Final verification of ROI filtering
clear all;
cd('/Users/jf5479/Dropbox/MATLAB/onPaths/YAHT');

fprintf('=== FINAL ROI FILTERING TEST ===\n\n');

% Test with all user's animals
theseAnimals = {'JF107', 'JF093', 'JF091', 'JF104', 'JF105', 'JF118', 'JF119', 'JF120', 'JF121'};
regionsNames = {'CP', 'GPe', 'SNr'};
plotHemisphere = [-1, 1, -1];

fprintf('Testing with %d animals\n', length(theseAnimals));
fprintf('Regions: %s\n', strjoin(regionsNames, ', '));
fprintf('Hemisphere settings: CP=left, GPe=right, SNr=left\n');
fprintf('onlyROIProbes = true (should show ONLY ROI portions)\n\n');

try
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, [], true, true, false, true, [], [], [], plotHemisphere);
    fprintf('\n✓ SUCCESS! Function completed with ROI filtering.\n');
    fprintf('The probes should now show ONLY the portions within CP, GPe, and SNr.\n');
    fprintf('NOT the full dorsal-ventral extent.\n');
catch e
    fprintf('\n✗ ERROR: %s\n', e.message);
    fprintf('At line %d\n', e.stack(1).line);
end