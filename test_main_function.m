% Test the main function with debugging
clear all;
cl_myPaths;

% Exact parameters from user's call
theseAnimals = {'JF107', 'JF093', 'JF091', 'JF104', 'JF105', 'JF118', 'JF119', 'JF120', 'JF121'};
regionsNames = {'CP', 'GPe', 'SNr'};
patchBrain = [];
onlyROIProbes = true;  % This is probably true based on the error
showPoints = true;
useBezierFit = false;
showRegionPlot = true;
regionColors = [];
blackBackground = [];
thickBrainLines = [];
plotHemisphere = [-1, 1, -1];

fprintf('Testing main function with onlyROIProbes = %d\n', onlyROIProbes);

% Call the function
try
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, patchBrain, onlyROIProbes, showPoints, ...
        useBezierFit, showRegionPlot, regionColors, blackBackground, thickBrainLines, plotHemisphere);
    fprintf('Function completed successfully\n');
catch e
    fprintf('Function failed: %s\n', e.message);
    disp(e.stack);
end