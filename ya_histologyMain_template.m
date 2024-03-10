
%% Histology main
clearvars -global % releases previous GUIs, if there are any

% EDIT ME :
animal = 'JF117';
brainsawPath = ['/home/netshare/zaru/', animal, '/Histology/downsampled_stacks/025_micron/'];
brainglobeLocation = '/home/julie/.brainglobe/'; % where your brainglobe data lives
allenAtlasPath = '/home/julie/Dropbox/Atlas/allenCCF'; % location of the allencff downloaded files

% registration parameters
orientationType = 'psl'; % psl (for posterior, superior, left), means the first, top left voxel
% is the most posterior, superior, left part of the brain
channelColToRegister = 'green'; % channel you want to register
channelColToTransform = 'red'; % channel you want to use to draw probes
atlasResolution_um = 25; % voxel size. currently in all dimensions, both for atlas and acquired image
atlasSpecies = 'mouse'; % atlas species
atlasType = 'allen'; % atlas name


%% ~ Images info ~
% registration location/files
[atlasLocation, imgToRegister, imgToTransform, outputDir] = ...
    ya_getLocations(brainglobeLocation, brainsawPath, channelColToRegister, ...
    channelColToTransform, atlasType, atlasSpecies, atlasResolution_um);

%% ~ Load in images and template ~
[tv, av, st, bregma] = ya_loadAllenAtlas([atlasLocation.folder, filesep, atlasLocation.name]);

%% ~ Register ~
if isempty(dir([outputDir, filesep, 'downsampled_standard.tiff']))
    % process image
    imgToRegister_processed = ya_preprocessImage([imgToRegister.folder, filesep, imgToRegister.name]);
    % register image
    ya_brainreg([imgToRegister_processed.folder, filesep, imgToRegister_processed.name], outputDir, ...
        [imgToTransform.folder, filesep, imgToTransform.name], orientationType, atlasResolution_um) % run brainreg
    % load image
    registeredImage = loadtiff([outputDir, filesep, 'downsampled_standard.tiff']);
    % generate appropriate formated files
    ya_convertToAPFormat(registeredImage, tv, av, outputDir) % get and save atlas data in standard
else
    % load image
    registeredImage = loadtiff([outputDir, filesep, 'downsampled_standard.tiff']);
end

%% ~ Manually check and adjust registration ~
%% ~~ Check and adjust orientation ~~
% Adjust the atlas to one slice (using arrow keys to navigate in the allen ccf and pressing eneter once you're happy)
% and this transformation will then be applied to all other slices and saved in a histology_ccf.mat file.
% If you register several slices, the average transform will be applied to
% all other slices and saved.
ya_checkAndCorrectAtlasOrientation(tv, av, st, registeredImage, outputDir);

%% ~~ Check adjust border alignement/scaling ~~
ya_checkAndCorrectAtlasAlignment(tv, av, st, registeredImage, outputDir)
% histology_ccf.mat : corresponding CCF slices
% atlas2histology_tform.mat : ccf to histology alignement/warp

%% ~ Draw probes + adjust probe fits ~
transformedImageDir = dir([outputDir, filesep, 'downsampled_standard_*.tiff']);
transformedImage = loadtiff([transformedImageDir.folder, filesep, transformedImageDir.name]);
ya_drawProbes_bezierCurves(tv, av, st, transformedImage, outputDir)

%% ~ Assign probes to days/sites ~
% For each probe, assign the day and recording site. If you have 4-shank
% probes also assign the shank number (from 1 to 4) - otherwise write NaN. 
% Optionally, add a 'certainty' value to keep track of any probes you are not
% 100% sure of. 

load([outputDir, '/probe_ccf.mat'])
ya_plotHistoPerMouse(outputDir, st);

% plot and assign
probe2ephys = struct;
% Probe 1 information
probe2ephys(1).day = 5;
probe2ephys(1).site = 1;
probe2ephys(1).shank = NaN;
probe2ephys(1).certainty = 0;
% Probe 2 information
probe2ephys(2).day = 4; %4
probe2ephys(2).site = 1;
probe2ephys(2).shank = NaN;
probe2ephys(2).certainty = 0;
% Probe 3 information
probe2ephys(3).day = 1; %1
probe2ephys(3).site = 1;
probe2ephys(3).shank = NaN;
probe2ephys(3).certainty = 0;
% Probe 4 information
probe2ephys(4).day = 3; %3
probe2ephys(4).site = 1;
probe2ephys(4).shank = NaN;
probe2ephys(4).certainty = 1;
% Probe n information
% .... and so on! 

% save 
save([outputDir, '/probe2ephys.mat'], 'probe2ephys')

%% ~ Align ephys and histology ~
load([outputDir, '/probe2ephys.mat']);
iProbe = 1;
%% Run this section for each probe 
% run this section for each probe, one at a time (otherwise the alignements
% will not be saved properly). Modify the loading function path input so it
% corresponds to how your data is stored 

iProbe = iProbe + 1;

% get information for this probe 
site = probe2ephys(iProbe).site;
day = probe2ephys(iProbe).day;
shank = probe2ephys(iProbe).shank;
pathToEphys = ['/home/netshare/zaru/', animal, filesep, num2str(day), filesep, 'ephys', filesep 'site']; % EDIT ME 
probeLength = 3840; % EDIT ME, size of recording sites in um. E.g: 3840 for NP1, 2880 for NP2, ect. 
ephys_sample_rate = 3000; % EDIT ME ephys sample rate (samples per second). 

% load data
[spike_times, spike_templates, template_depths, spike_depths, template_xdepths,...
    spike_xdepths, channelPositions, goodChannels] = ya_loadEphysData(pathToEphys, ephys_sample_rate);

% align ephys and histology. Drag the lines to stretch or shrink any
% regions that need it
ya_alignEphysAndHistology_draggable(st, outputDir, ...
    spike_times, spike_templates, template_depths, ...
    spike_xdepths, template_xdepths, lfp, channel_positions(:, 2), channel_positions(:, 1), ...
    iProbe, probeLength, shank);

%% ~ Find out which units are in which region / which anatomical position ~
load([outputDir, '/probe2ephys.mat'])
load([outputDir, '/probe_ccf.mat'])

iProbe = 1;
iUnit = 1;
% get 
unit_idx = find(probe_ccf(iProbe).probe_depths >= template_depths(iUnit), 1, 'first');
% ML * DV * AP coordinates 
unit_coords = probe_ccf(iProbe).trajectory_coords(unit_idx,:);

% area 
unit_region = st.acronym(find(st.id == probe_ccf(iProbe).trajectory_areas(unit_idx,:)));

%% ~ Various useful plotting functions ~
% plot the probe tracks in the current mouse
ya_plotHistoPerMouse(outputDir, st)

% plot/generate a video of all your probe tracks across mice
paths = {[outputDir, '/probe_ccf.mat'], [outputDir, '/probe_ccf.mat']};
ya_plotAllProbeTracksInROIs(allenAtlasPath, [atlasLocation.folder, filesep, atlasLocation.name], paths, regionNames, '', '') 

% plot/generate a video all probe tracks in certain regions, across mice
regionNames = {'CP', 'GPe', 'SNr'};
regionColors = [0, 0.7461, 1; 0.1797, 0.5430, 0.3398;1, 0.4102, 0.7031];
regionPlotSide = [-1, 1, -1];
ya_plotAllProbeTracksInROIs(allenAtlasPath, [atlasLocation.folder, filesep, atlasLocation.name],...
    paths, regionNames, regionColors, regionPlotSide) 

% plot a brain outline with certain brain regions
ya_plotAllProbeTracksInROIs(allenAtlasPath, [atlasLocation.folder, filesep, atlasLocation.name],...
    {}, regionNames, regionColors, regionPlotSide);

% use a patch outline rather than grid
plotPatch = 1;
ya_plotAllProbeTracksInROIs(allenAtlasPath, [atlasLocation.folder, filesep, atlasLocation.name],...
    {}, regionNames, regionColors, regionPlotSide, plotPatch);

