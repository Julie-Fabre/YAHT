% Diagnose key variables that the function needs
clear all;

fprintf('=== Checking required variables and paths ===\n');

% Run cl_myPaths
try
    cl_myPaths;
    fprintf('✓ cl_myPaths executed successfully\n');
catch e
    fprintf('✗ cl_myPaths failed: %s\n', e.message);
    return;
end

% Check if key variables are defined
fprintf('\n=== Checking key variables ===\n');

if exist('allenAtlasPath', 'var')
    fprintf('✓ allenAtlasPath = %s\n', allenAtlasPath);
    if exist(allenAtlasPath, 'dir')
        fprintf('  ✓ Directory exists\n');
    else
        fprintf('  ✗ Directory does not exist!\n');
    end
else
    fprintf('✗ allenAtlasPath not defined!\n');
end

if exist('brainglobeLocation', 'var')
    fprintf('✓ brainglobeLocation = %s\n', brainglobeLocation);
    if exist(brainglobeLocation, 'dir')
        fprintf('  ✓ Directory exists\n');
    else
        fprintf('  ✗ Directory does not exist!\n');
    end
else
    fprintf('✗ brainglobeLocation not defined!\n');
end

% Check specific files that the function needs
fprintf('\n=== Checking required files ===\n');

if exist('allenAtlasPath', 'var')
    files_to_check = {
        'template_volume_10um.npy';
        'annotation_volume_10um_by_index.npy';
        'structure_tree_safe_2017.csv'
    };
    
    for i = 1:length(files_to_check)
        filepath = [allenAtlasPath, filesep, files_to_check{i}];
        if exist(filepath, 'file')
            fprintf('  ✓ %s\n', files_to_check{i});
        else
            fprintf('  ✗ %s MISSING!\n', files_to_check{i});
        end
    end
end

% Check brainglobe atlas
fprintf('\n=== Checking brainglobe atlas ===\n');
if exist('brainglobeLocation', 'var')
    atlasType = 'allen';
    atlasSpecies = 'mouse';
    atlasResolution_um = 25;
    atlasLocation = dir([brainglobeLocation, atlasType, '_', atlasSpecies, '_', num2str(atlasResolution_um), 'um*']);
    
    if ~isempty(atlasLocation)
        fprintf('  ✓ Found brainglobe atlas: %s\n', atlasLocation.name);
    else
        fprintf('  ✗ Brainglobe atlas not found!\n');
    end
end

% Test a sample animal file
fprintf('\n=== Checking sample animal file ===\n');
testAnimal = 'JF107';
outputDir = ['/Users/jf5479/Dropbox/Histology/' testAnimal '/brainReg/'];
probe_file = [outputDir, 'probe_ccf.mat'];

if exist(probe_file, 'file')
    fprintf('  ✓ %s probe file exists\n', testAnimal);
    try
        load(probe_file);
        fprintf('  ✓ Loaded %d probes\n', length(probe_ccf));
        
        % Check first probe structure
        if length(probe_ccf) > 0
            hasPoints = isfield(probe_ccf(1), 'points') && ~isempty(probe_ccf(1).points);
            hasTrajectory = isfield(probe_ccf(1), 'trajectory_areas') && ~isempty(probe_ccf(1).trajectory_areas);
            fprintf('  ✓ First probe: points=%d, trajectory=%d\n', hasPoints, hasTrajectory);
        end
    catch e
        fprintf('  ✗ Failed to load probe file: %s\n', e.message);
    end
else
    fprintf('  ✗ %s probe file missing!\n', testAnimal);
end

fprintf('\n=== Diagnosis complete ===\n');