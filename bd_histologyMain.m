%% ~~ Histology main ~~ %% 

%% Images info %% 
myPaths; % see JF_scripts_cortexlab
animal = 'JF070';
orientationType = 'psl';
channelColToRegister = 'green';
atlas = 25; 
atlasLocation = '/home/julie/.brainglobe/allen_mouse_25um_v1.2';

%% Load in images and template %% 
%bd_loadAllenAtlas()
% load in allen atlas: get in 'AP' format: histology_ccf.mat with tv_slices, av_slices,
% plane_ap, plane_ml, plane_dv 
allen_atlas_path = [allenAtlasPath 'allenCCF'];
% tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
% av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
tv = loadtiff([atlasLocation, filesep, 'reference.tiff']);
av = loadtiff([atlasLocation, filesep, 'annotation.tiff']);
st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
bregma = [540,0,570];
%allenAtlas10um = readNPY([allenAtlasPath 'allenCCF' filesep 'template_volume_10um.npy']);

% load in channel to register 
channelToRegister = dir([brainsawPath, '/*/', animal, '/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
outputDir = [channelToRegister.folder, filesep, 'brainReg'];

%% Register %% 
bd_brainreg([channelToRegister.folder, filesep, channelToRegister.name], outputDir, orientationType , atlas)
registeredImage = loadtiff([outputDir, filesep, 'downsampled_standard.tiff']);
bd_convertToAPFormat(registeredImage, tv, av, outputDir)

%% Manually check and adjust registration %% 
bd_checkAndCorrectAlign(tv,av,st,registeredImage,outputDir)

%% Draw probes %% 

%% Save data %%