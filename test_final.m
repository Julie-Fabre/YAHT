% Final test of the corrected function
clear all;
cl_myPaths;

% Test parameters
theseAnimals = {'JF107'};  % Start with one animal
regionsNames = {'CP', 'GPe', 'SNr'};
onlyROIProbes = true;  % Test ROI filtering
plotHemisphere = [-1, 1, -1];

fprintf('Testing corrected function with:\n');
fprintf('  Animals: %s\n', strjoin(theseAnimals, ', '));
fprintf('  Regions: %s\n', strjoin(regionsNames, ', '));
fprintf('  onlyROIProbes: %d\n', onlyROIProbes);
fprintf('  plotHemisphere: [%s]\n', sprintf('%d ', plotHemisphere));

try
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, [], onlyROIProbes, true, ...
        false, false, [], [], [], plotHemisphere);
    fprintf('\n=== SUCCESS: Function completed without errors! ===\n');
catch e
    fprintf('\n=== ERROR: %s ===\n', e.message);
    if ~isempty(e.stack)
        fprintf('At line %d in %s\n', e.stack(1).line, e.stack(1).name);
    end
end