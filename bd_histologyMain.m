
%% Histology main
% - add option to use elastix
% - don't fit some parts that are prone to being moved (eg olf bulbs, pposterior cortex (eg retrospenial
% ect) ?

%% ~ Images info
myPaths; % see JF_scripts_cortexlab
animal = 'JF070';

% registration parameters
orientationType = 'psl';
channelColToRegister = 'green';
channelColToTransform = 'red';
atlasResolution_um = 25;
atlasSpecies = 'mouse';
atlasType = 'allen';
atlasSize = [320, 456, 528]; 

% registration location/files
atlasLocation = dir(['/home/julie/.brainglobe/', atlasType, '_', atlasSpecies, '_', num2str(atlasResolution_um), 'um*']); % atlas location + atlas to use
imgToRegister = dir([brainsawPath, '/*/', animal, '/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
imgToTransform = dir([brainsawPath, '/*/', animal, '/downsampled_stacks/025_micron/*', channelColToTransform, '*.tif*']);
outputDir = [imgToRegister.folder, filesep, 'brainReg'];

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

%% ~~ Check adjust border alignement/scaling ~~
bd_checkAndCorrectAlign(tv, av, st, registeredImage, outputDir, screenToUse)
% histology_ccf.mat : corresponding CCF slices
% atlas2histology_tform.mat : histology/ccf alignement/warp

%% ~ Draw probes ~
bd_drawProbes(tv, av, st, registeredImage, outputDir)

%% ~ Assign probes to days/sites and save data ~