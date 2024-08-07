
%% Histology main
% To do:
% - add option to use elastix
% - don't fit some parts that are prone to being moved (eg olf bulbs, posterior cortex (eg retrospenial
% ect) ?
% - add option to change atlas orientation post automatic alignement
% - two views in probe drawing GUI
% - GUI to assign probes

%% ~ Images info
clear all;
close all;

clearvars -global % releases previous GUIs, if there are any
cl_myPaths; % see https://github.com/Julie-Fabre/JF_Scripts_CortexLab/blob/master/load/myPaths.m.
% loads in a bunch of paths, of which only one is used in this script: brainsawPath
animal = 'JF109';

% registration parameters
orientationType = 'psl'; % psl (for posterior, superior, left), means the first, top left voxel
% is the most posterior, superior, left part of the brain
channelColToRegister = 'green'; % channel you want to register
channelColToTransform = 'red'; % channel you want to use to draw probes
atlasResolution_um = 25; % voxel size. currently in all dimensions, both for atlas and acquired image
atlasSpecies = 'mouse'; % atlas species
atlasType = 'allen'; % atlas name
brainglobeLocation = '/home/julie/.brainglobe/'; % where your brainglobe data lives

brainsawPath_curr = [cl_cortexlab_filename(animal, '', '', 'histo_folder', '', '', ''), '/downsampled_stacks/025_micron'];
% registration location/files
[atlasLocation, imgToRegister, imgToTransform, outputDir] = ...
    ya_getLocations(brainglobeLocation, brainsawPath_curr, channelColToRegister, ...
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
screenToUse = 2; % on which of your displays to create the following plots. 1 = main display, 2 = secondary and so on

%% ~~ Check and adjust orientation ~~
% Adjust the atlas to one slice (using arrow keys to navigate in the allen ccf and pressing eneter once you're happy)
% and this transformation will then be applied to all other slices and saved in a histology_ccf.mat file.
% If you register several slices, the average transform will be applied to
% all other slices and saved.
if isempty(dir([outputDir, filesep, 'manual', filesep, 'histology_ccf.mat']))
    ya_checkAndCorrectAtlasOrientation(tv, av, st, registeredImage, outputDir, screenToUse);
end

%% ~~ Check adjust border alignement/scaling ~~
if isempty(dir([outputDir, filesep, 'manual', filesep, 'atlas2histology_tform.mat']))
    ya_checkAndCorrectAtlasAlignment(tv, av, st, registeredImage, outputDir, screenToUse)
end
% histology_ccf.mat : corresponding CCF slices
% atlas2histology_tform.mat : ccf to histology alignement/warp

%% ~ Draw probes + adjust probe fits ~
transformedImageDir = dir([outputDir, filesep, 'downsampled_standard_*.tiff']);
transformedImage = loadtiff([transformedImageDir.folder, filesep, transformedImageDir.name]);
screenOri = 'portrait';
ya_drawProbes_bezierCurves(tv, av, st, transformedImage, outputDir, screenOri)

load([outputDir, '/probe_ccf.mat'])

save([outputDir, '/probe_ccf.mat'], 'probe_ccf')
save([outputDir, '/probe_ccf', date, '1.mat'], 'probe_ccf')

%% ~ Assign probes to days/sites ~
load([outputDir, '/probe_ccf.mat'])
ya_plotHistoPerMouse(outputDir, st);

% duplicate probes if several recordings (eg botRow+0, botRow+48, ...)
probe_ccf(end+1) = probe_ccf(1);
probe_ccf(end+1) = probe_ccf(2);
probe_ccf(end+1) = probe_ccf(3);
probe_ccf(end+1) = probe_ccf(4);

probe_ccf(end+1) = probe_ccf(1);
probe_ccf(end+1) = probe_ccf(2);
probe_ccf(end+1) = probe_ccf(3);
probe_ccf(end+1) = probe_ccf(4);

probe_ccf(end+1) = probe_ccf(7);
probe_ccf(end+1) = probe_ccf(8);
probe_ccf(end+1) = probe_ccf(9);
probe_ccf(end+1) = probe_ccf(10);

probe_ccf(end+1) = probe_ccf(7);
probe_ccf(end+1) = probe_ccf(8);
probe_ccf(end+1) = probe_ccf(9);
probe_ccf(end+1) = probe_ccf(10);

% plot and assign
probe2ephys = struct;
probe2ephys(1).day = 4;
probe2ephys(1).site = 1;
probe2ephys(1).shank = NaN;
probe2ephys(1).certainty = 0;

probe2ephys(2).day = 6; %4
probe2ephys(2).site = 1;
probe2ephys(2).shank = NaN;
probe2ephys(2).certainty = 0;

probe2ephys(3).day = 1; %1
probe2ephys(3).site = 1;
probe2ephys(3).shank = 4;
probe2ephys(3).certainty = 0;

probe2ephys(4).day = 1; %3
probe2ephys(4).site = 1;
probe2ephys(4).shank = 3;
probe2ephys(4).certainty = 1;

probe2ephys(5).day = 5;
probe2ephys(5).site = 1;
probe2ephys(5).shank = NaN;
probe2ephys(5).certainty = 0;

probe2ephys(6).day = 2;
probe2ephys(6).site = 1;
probe2ephys(6).shank = 4;
probe2ephys(6).certainty = 1;

probe2ephys(7).day = 7;
probe2ephys(7).site = 1;
probe2ephys(7).shank = NaN;
probe2ephys(7).certainty = 1;


probe2ephys(8).day = 3;
probe2ephys(8).site = 1;
probe2ephys(8).shank = 2;
probe2ephys(8).certainty = 1;

probe2ephys(9).day = 3;
probe2ephys(9).site = 1;
probe2ephys(9).shank = 1;
probe2ephys(9).certainty = 1;

probe2ephys(10).day = 2;
probe2ephys(10).site = 1;
probe2ephys(10).shank = 1;
probe2ephys(10).certainty = 1;

probe2ephys(11).day = 1;
probe2ephys(11).site = 1;
probe2ephys(11).shank = 1;
probe2ephys(11).certainty = 1;
%
probe2ephys(12).day = 3;
probe2ephys(12).site = 1;
probe2ephys(12).shank = 3;
probe2ephys(12).certainty = 0;

probe2ephys(13).day = 2;
probe2ephys(13).site = 1;
probe2ephys(13).shank = 2;
probe2ephys(13).certainty = 0;

probe2ephys(14).day = 2;
probe2ephys(14).site = 1;
probe2ephys(14).shank = 3;
probe2ephys(14).certainty = 0;

probe2ephys(15).day = 3;
probe2ephys(15).site = 1;
probe2ephys(15).shank = 4;
probe2ephys(15).certainty = 0;

probe2ephys(16).day = 1;
probe2ephys(16).site = 1;
probe2ephys(16).shank = 2;
probe2ephys(16).certainty = 0;

probe2ephys(17).day = 1;
probe2ephys(17).site = 2;
probe2ephys(17).shank = NaN;
probe2ephys(17).certainty = 0;

probe2ephys(18).day = 1;
probe2ephys(18).site = 3;
probe2ephys(18).shank = NaN;
probe2ephys(18).certainty = 0;

probe2ephys(19).day = 7;
probe2ephys(19).site = 2;
probe2ephys(19).shank = NaN;
probe2ephys(19).certainty = 0;

probe2ephys(20).day = 2;
probe2ephys(20).site = 2;
probe2ephys(20).shank = NaN;
probe2ephys(20).certainty = 0;

probe2ephys(21).day = 3;
probe2ephys(21).site = 2;
probe2ephys(21).shank = NaN;
probe2ephys(21).certainty = 0;

probe2ephys(22).day = 2;
probe2ephys(22).site = 3;
probe2ephys(22).shank = NaN;
probe2ephys(22).certainty = 0;

probe2ephys(23).day = 3;
probe2ephys(23).site = 3;
probe2ephys(23).shank = NaN;
probe2ephys(23).certainty = 0;
%
probe2ephys(24).day = 5;
probe2ephys(24).site = 2;
probe2ephys(24).shank = NaN;
probe2ephys(24).certainty = 0;

probe2ephys(25).day = 4;
probe2ephys(25).site = 2;
probe2ephys(25).shank = NaN;
probe2ephys(25).certainty = 0;

probe2ephys(26).day = 6;
probe2ephys(26).site = 2;
probe2ephys(26).shank = NaN;
probe2ephys(26).certainty = 0;

probe2ephys(27).day = 8;
probe2ephys(27).site = 1;
probe2ephys(27).shank = NaN;
probe2ephys(27).certainty = 0;

probe2ephys(28).day = 6;
probe2ephys(28).site = 2;
probe2ephys(28).shank = NaN;
probe2ephys(28).certainty = 0;

save([outputDir, '/probe2ephys.mat'], 'probe2ephys')
save([outputDir, '/probe2ephys', date, '.mat'], 'probe2ephys')

save([outputDir, '/probe_ccf.mat'], 'probe_ccf')
save([outputDir, '/probe_ccf', date, '.mat'], 'probe_ccf')

%%
probe_points2 = probe_points;
nonEmptyCells = find(~cellfun(@isempty, probe_points2(:, 16)));

for iNonEmptyCells = 1:size(nonEmptyCells, 1)
    thisCell = nonEmptyCells(iNonEmptyCells);
    probe_points2(thisCell-18, 13) = probe_points2(thisCell, 16);
end


probe_points = probe_points2;
save([outputDir, '/probe_points.mat'], 'probe_points')

%% ~ Align ephys and histology ~
% - add "5 min histology recs" look at
% - check out imro file , plot
% - plot probe in image overlay (rotated on probe axis)
iProbe = 0;
load([outputDir, '/probe2ephys.mat'])
load([outputDir, '/probe_ccf.mat'])

%%
%for iProbe = 1:17
%$for iProbe = 5%[20,21,22,23,26,27]
for iProbe = 6%[2,3,4,10,11]

    keep probe2ephys animal iProbe st outputDir tv av

    % experiment info
    site = probe2ephys(iProbe).site; %if ~isnan(site)
    experiments = cl_find_experiments(animal, '', true);
    experiments = experiments([experiments.ephys]);
    day_num = probe2ephys(iProbe).day;
    thisDate = experiments(day_num).thisDate;
    shank = probe2ephys(iProbe).shank;

    % select experiment number
    %if site == 1 || site == 2
    experiment = 3;
    %elseif site == 3 || site == 4
    %    experiment = 1;
    %elseif site == 5 || site == 6
    %    experiment = 5;
    %end


    % recording
    if isfield(probe2ephys, 'recording') && ~isempty(probe2ephys(iProbe).recording)
        recording = probe2ephys(iProbe).recording(1);
        experiment = 5;
    else
        recording = [];

    end

    try
        experiment= 2
        % load ephys data
        load_sync = true;
        load_parts.cam = false;
        loadClusters = 0;
        load_parts.ephys = true;
        cl_load_experiment;

        % plot PSTH
        curr_shank = shank;
        lfp = NaN;
        cl_cellraster({stimOn_times}, {trial_conditions(:, 1)})
    catch
        try
            experiment = 3;
            % load ephys data
            load_sync = true;
            load_parts.cam = false;
            loadClusters = 0;
            load_parts.ephys = true;
            cl_load_experiment;

            % plot PSTH
            curr_shank = shank;
            lfp = NaN;
            cl_cellraster({stimOn_times}, {trial_conditions(:, 1)})
        catch
            experiment = 1;
            % load ephys data
            load_sync = true;
            load_parts.cam = false;
            loadClusters = 0;
            load_parts.ephys = true;
            cl_load_experiment;

            % plot PSTH
            curr_shank = shank;
            lfp = NaN;
            cl_cellraster({stimOn_times}, {trial_conditions(:, 1)})
        end
    end
    % get probe length specs
    if max(template_depths) > 2880
        probeLength = 3840;
    elseif max(template_depths) > 1500
        probeLength = 2880;
    elseif max(template_depths) > 750
        probeLength = 1500;
    elseif max(template_depths) <= 750
        probeLength = 750;
    end
    try
        trial_conditions_clean = trial_conditions;
        trial_conditions_clean(trial_conditions(:, 1) == 4, 1) = 100; %go1
        trial_conditions_clean(trial_conditions(:, 1) == 12, 1) = 101; %go2
        trial_conditions_clean(trial_conditions(:, 1) == 6, 1) = 102; %no go
        %trial_conditions_clean(trial_conditions(:,1)==11,1) = 103;%no like
        %trial_conditions_clean(trial_conditions(:,1)==13,1) = 104;%go like
        theseImages = [100:102];


        theseImages_trials = ismember(trial_conditions_clean(:, 1), theseImages) & ismember(trial_conditions_clean(:, 2), [-90]);
        cl_cellraster({stimOn_times(theseImages_trials), stimOn_times(theseImages_trials), stimOn_times(theseImages_trials)}, ...
            {trial_conditions_clean(theseImages_trials, 1), ...
            trial_conditions_clean(theseImages_trials, 2), -89 + trial_conditions_clean(theseImages_trials, 1) + abs(trial_conditions_clean(theseImages_trials, 2))})
    end

    % get quality metrics
    try
        rerunQM = 0;

        [unitType, qMetric] = bc_qualityMetricsPipeline_JF(animal, thisDate, site, recording, experiment, [], rerunQM, 0, 1);
    catch
        unitType = zeros(size(template_depths, 1), 1);
    end
    % align
    ya_alignEphysAndHistology_draggable(st, outputDir, ...
        spike_times, spike_templates, template_depths, ...
        spike_xdepths, template_xdepths, lfp, channel_positions(:, 2), channel_positions(:, 1), ...
        iProbe, probeLength, shank, unitType);
    %end
    % copy to backup
    % copyfile([outputDir, '/probe_ccf.mat'],['/media/julie/Expansion/Histology_backup/', animal, '/probe_ccf.mat'])
end

%% final
copyfile([outputDir, '/probe_ccf.mat'], ['/media/julie/ExtraHD/Histology_backup/', animal, '/probe_ccf.mat'])

copyfile([outputDir, '/probe2ephys.mat'], ['/media/julie/ExtraHD/Histology_backup/', animal, '/probe2ephys.mat'])

%%
load([outputDir, '/probe_ccf.mat'])
probe_ccf_new = probe_ccf;
load(['/media/julie/Expansion/Histology_backup/', animal, '/probe_ccf.mat'])
probe_ccf(1) = probe_ccf_new(1);
probe_ccf(2) = probe_ccf_new(2);
probe_ccf(5) = probe_ccf_new(5);
save([outputDir, '/probe_ccf.mat'], 'probe_ccf')

%% ~ Various useful plotting functions ~