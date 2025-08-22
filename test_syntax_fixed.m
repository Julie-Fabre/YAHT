% Test after fixing syntax error
clear all;
cl_myPaths;

fprintf('Testing function after syntax fix...\n');

% Test with minimal parameters first
theseAnimals = {'JF107'};
regionsNames = {'CP', 'GPe', 'SNr'};
onlyROIProbes = false;  % Start with false to see all probes
plotHemisphere = [-1, 1, -1];

try
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, [], onlyROIProbes, true, ...
        false, false, [], [], [], plotHemisphere);
    fprintf('SUCCESS: Function ran without syntax errors!\n');
catch e
    fprintf('ERROR: %s\n', e.message);
    if ~isempty(e.stack)
        fprintf('Line %d in %s\n', e.stack(1).line, e.stack(1).name);
    end
end