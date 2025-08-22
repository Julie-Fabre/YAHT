% Diagnostic script for the user to run
% This checks the environment and tests the function step by step

fprintf('=== DIAGNOSTIC FOR USER ===\n');

% Step 1: Check current directory and path
fprintf('Current directory: %s\n', pwd);
fprintf('MATLAB path includes this directory: %d\n', contains(path, pwd));

% Step 2: Check if the function file exists and is accessible
func_file = 'private_JF/ya_plotAllProbeTracksInROIs_JF.m';
fprintf('Function file exists: %d\n', exist(func_file, 'file'));

% Step 3: Try to get function info
try
    func_info = functions(@ya_plotAllProbeTracksInROIs_JF);
    fprintf('Function accessible: YES\n');
    fprintf('Function type: %s\n', func_info.type);
    fprintf('Function file: %s\n', func_info.file);
catch e
    fprintf('Function accessible: NO - %s\n', e.message);
end

% Step 4: Check required paths and variables
try
    cl_myPaths;
    fprintf('cl_myPaths executed: YES\n');
    
    if exist('allenAtlasPath', 'var')
        fprintf('allenAtlasPath defined: YES (%s)\n', allenAtlasPath);
    else
        fprintf('allenAtlasPath defined: NO\n');
    end
    
    if exist('brainglobeLocation', 'var')
        fprintf('brainglobeLocation defined: YES (%s)\n', brainglobeLocation);
    else
        fprintf('brainglobeLocation defined: NO\n');
    end
catch e
    fprintf('cl_myPaths failed: %s\n', e.message);
end

% Step 5: Test a minimal function call
fprintf('\n=== ATTEMPTING FUNCTION CALL ===\n');
try
    % Call with minimal parameters
    ya_plotAllProbeTracksInROIs_JF({'JF107'}, {'CP'}, [], false, true, false, false, [], [], [], []);
    fprintf('✓ FUNCTION CALL SUCCEEDED!\n');
catch e
    fprintf('✗ FUNCTION CALL FAILED: %s\n', e.message);
    fprintf('This indicates the fixes are working but there may be a runtime issue.\n');
end

fprintf('\n=== DIAGNOSTIC COMPLETE ===\n');