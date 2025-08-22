% Test with user's exact parameters
clear all;
cl_myPaths;

fprintf('Testing with user parameters...\n');

% User's exact call parameters
theseAnimals = {'JF107', 'JF093', 'JF091', 'JF104', 'JF105', 'JF118', 'JF119', 'JF120', 'JF121'};
regionsNames = {'CP', 'GPe', 'SNr'};
onlyROIProbes = true;  % This was causing the "no probes" issue
showPoints = true;
useBezierFit = false;
showRegionPlot = true;
regionColors = [];
blackBackground = [];
thickBrainLines = [];
plotHemisphere = [-1, 1, -1];

fprintf('Animals: %s\n', strjoin(theseAnimals, ', '));
fprintf('Regions: %s\n', strjoin(regionsNames, ', '));
fprintf('onlyROIProbes: %d\n', onlyROIProbes);

try
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, [], onlyROIProbes, showPoints, ...
        useBezierFit, showRegionPlot, regionColors, blackBackground, thickBrainLines, plotHemisphere);
    fprintf('SUCCESS: Function completed with user parameters!\n');
catch e
    fprintf('ERROR: %s\n', e.message);
    if ~isempty(e.stack)
        fprintf('Line %d in %s\n', e.stack(1).line, e.stack(1).name);
    end
end