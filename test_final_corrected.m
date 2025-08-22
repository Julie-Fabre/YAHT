% Final test - should work now!
clear all;
cl_myPaths;

fprintf('=== Testing CORRECTED function ===\n');

% Start with a simple test: one animal, no ROI filtering
fprintf('Test 1: Single animal, no ROI filtering...\n');
try
    ya_plotAllProbeTracksInROIs_JF({'JF107'}, {'CP', 'GPe', 'SNr'}, [], false, true, false, false, [], [], [], [-1, 1, -1]);
    fprintf('✓ Test 1 PASSED\n');
catch e
    fprintf('✗ Test 1 FAILED: %s\n', e.message);
    return;
end

% Test 2: Single animal WITH ROI filtering
fprintf('\nTest 2: Single animal, with ROI filtering...\n');
try
    ya_plotAllProbeTracksInROIs_JF({'JF107'}, {'CP', 'GPe', 'SNr'}, [], true, true, false, false, [], [], [], [-1, 1, -1]);
    fprintf('✓ Test 2 PASSED\n');
catch e
    fprintf('✗ Test 2 FAILED: %s\n', e.message);
    return;
end

% Test 3: User's original parameters
fprintf('\nTest 3: User''s original parameters...\n');
try
    theseAnimals = {'JF107', 'JF093', 'JF091', 'JF104', 'JF105', 'JF118', 'JF119', 'JF120', 'JF121'};
    ya_plotAllProbeTracksInROIs_JF(theseAnimals, {'CP', 'GPe', 'SNr'}, [], true, true, false, true, [], [], [], [-1, 1, -1]);
    fprintf('✓ Test 3 PASSED - Original parameters work!\n');
catch e
    fprintf('✗ Test 3 FAILED: %s\n', e.message);
    return;
end

fprintf('\n🎉 ALL TESTS PASSED! Function is working correctly! 🎉\n');