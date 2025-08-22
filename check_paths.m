% Check what cl_cortexlab_filename returns vs hardcoded path
clear all;
cl_myPaths;

animal = 'JF107';

% Method 1: cl_cortexlab_filename (used by plotProbeRegions)
path1 = cl_cortexlab_filename(animal, [], [], 'histo');
fprintf('cl_cortexlab_filename: %s\n', path1);
fprintf('  exists: %d\n', exist(path1, 'file'));

% Method 2: hardcoded path (used by main function)
outputDir = ['/Users/jf5479/Dropbox/Histology/' animal '/brainReg/'];
path2 = [outputDir, 'probe_ccf.mat'];
fprintf('hardcoded path: %s\n', path2);
fprintf('  exists: %d\n', exist(path2, 'file'));