%% ~~ Histology main ~~ %% 

%% Images info %% 
myPaths; % see JF_scripts_cortexlab
animal = 'JF070';

%% Load in images and template %% 
% get in 'AP' format: histology_ccf.mat with tv_slices, av_slices,
% plane_ap, plane_ml, plane_dv 
allen_atlas_path = [allenAtlasPath 'allenCCF'];
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
bregma = [540,0,570];
allenAtlas10um = readNPY([allenAtlasPath 'allenCCF' filesep 'template_volume_10um.npy']);

channelColToRegister = 'green';
channelToRegister = dir([brainsawPath, '/*/', animal, '/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
orientationType = 'psl';
outputDir = channelToRegister.folder;

%% Register %% 
bd_brainreg([channelToRegister.folder, filesep, channelToRegister.name], outputDir, 'psl')
%% Manually check and adjust registration %% 

%% Draw probes %% 

%% Save data %%