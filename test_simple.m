% Simple test to trace execution
clear all;
cl_myPaths;

% Test with minimal parameters
theseAnimals = {'JF107'};  % Just one animal to start
regionsNames = {'CP', 'GPe', 'SNr'};
onlyROIProbes = false;  % Start with false to see if any probes are detected
showPoints = true;
useBezierFit = false;
showRegionPlot = false;  % Disable region plot to focus on main plotting
plotHemisphere = [-1, 1, -1];

fprintf('=== Testing with single animal and onlyROIProbes=false ===\n');

try
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, [], onlyROIProbes, showPoints, ...
        useBezierFit, showRegionPlot, [], [], [], plotHemisphere);
    fprintf('SUCCESS: Function completed\n');
catch e
    fprintf('ERROR: %s\n', e.message);
    if ~isempty(e.stack)
        fprintf('At line %d in %s\n', e.stack(1).line, e.stack(1).name);
    end
end