
%% Histology main
% - add option to use elastix
% - don't fit some parts that are prone to being moved (eg olf bulbs, pposterior cortex (eg retrospenial
% ect) ?

%% ~ Images info
myPaths; % see JF_scripts_cortexlab. loads in a bunch of paths, of which only one is used in this script: brainsawPath
animal = 'JF070';

% registration parameters
orientationType = 'psl'; % psl 
channelColToRegister = 'green'; % channel you want to register 
channelColToTransform = 'red'; % channel you want to use to draw probes
atlasResolution_um = 25; % voxel size. currently in all dimensions, both for atlas and acquired image
atlasSpecies = 'mouse'; % atlas species
atlasType = 'allen'; % atlas name
brainglobeLocation = '/home/julie/.brainglobe/'; % where your brainglobe data lives 

% registration location/files
[atlasLocation, imgToRegister, imgToTransform, outputDir] =...
    bd_getLocations(brainglobeLocation, atlasType, atlasSpecies, atlasResolution_um);

%% ~ Load in images and template ~
[tv, av, st, bregma] = bd_loadAllenAtlas([atlasLocation.folder, filesep, atlasLocation.name]);

%% ~ Register ~
bd_brainreg([imgToRegister.folder, filesep, imgToRegister.name], outputDir, ...
    [imgToTransform.folder, filesep, imgToTransform.name],orientationType, atlasResolution_um)
registeredImage = loadtiff([outputDir, filesep, 'downsampled_standard.tiff']);
bd_convertToAPFormat(registeredImage, tv, av, outputDir)

%% ~ Manually check and adjust registration ~
screenToUse = 2;

%% ~~ Check and adjust orientation ~~
bd_checkAndCorrectOrientation(tv, av, st, registeredImage, screenToUse);

%% ~~ [WIP] Check adjust border alignement/scaling ~~
bd_checkAndCorrectAlign(tv, av, st, registeredImage, outputDir, screenToUse)
% histology_ccf.mat : corresponding CCF slices
% atlas2histology_tform.mat : ccf to histology alignement/warp

%% ~ Draw probes ~
% QQ to add: draw both bspine fit and affine, choose which to use for each
% probe 
transformedImageDir = dir([outputDir, filesep, 'downsampled_standard_*.tiff']);
transformedImage = loadtiff([transformedImageDir.folder, filesep, transformedImageDir.name]);
bd_drawProbes(tv, av, st, transformedImage, outputDir)

%% ~ [WIP] Assign probes to days/sites and save data ~