% Test ROI filtering implementation
clear all;
cd('/Users/jf5479/Dropbox/MATLAB/onPaths/YAHT');

fprintf('=== Testing ROI Filtering ===\n');

% Test with one animal first
theseAnimals = {'JF107'};
regionsNames = {'CP', 'GPe', 'SNr'};
plotHemisphere = [-1, 1, -1];

% Test 1: Without ROI filtering (should show full probes)
fprintf('\nTest 1: onlyROIProbes = false (full probes)...\n');
try
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, [], false, true, false, false, [], [], [], plotHemisphere);
    fprintf('✓ Test 1 completed\n');
    pause(2);
    close all;
catch e
    fprintf('✗ Test 1 failed: %s\n', e.message);
    fprintf('  Line %d in %s\n', e.stack(1).line, e.stack(1).name);
end

% Test 2: With ROI filtering (should show only ROI portions)
fprintf('\nTest 2: onlyROIProbes = true (only ROI portions)...\n');
try
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, regionsNames, [], true, true, false, false, [], [], [], plotHemisphere);
    fprintf('✓ Test 2 completed\n');
    pause(2);
    close all;
catch e
    fprintf('✗ Test 2 failed: %s\n', e.message);
    fprintf('  Line %d in %s\n', e.stack(1).line, e.stack(1).name);
end

% Test 3: Full user parameters
fprintf('\nTest 3: All animals with ROI filtering...\n');
allAnimals = {'JF107', 'JF093', 'JF091', 'JF104', 'JF105', 'JF118', 'JF119', 'JF120', 'JF121'};
try
    ya_plotAllProbeTracksInROIs_JF(allAnimals, regionsNames, [], true, true, false, true, [], [], [], plotHemisphere);
    fprintf('✓ Test 3 completed - ROI filtering should be working!\n');
catch e
    fprintf('✗ Test 3 failed: %s\n', e.message);
    fprintf('  Line %d in %s\n', e.stack(1).line, e.stack(1).name);
end

fprintf('\n=== Testing Complete ===\n');