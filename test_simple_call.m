% Simple test that just calls the function
clear all;
cd('/Users/jf5479/Dropbox/MATLAB/onPaths/YAHT');
addpath('private_JF');

fprintf('Testing function call...\n');

% Set up variables first
theseAnimals = {'JF107'};
regionsNames = {'CP', 'GPe', 'SNr'};
plotHemisphere = [-1, 1, -1];

% Try the function call
ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, [], false, true, false, false, [], [], [], plotHemisphere);