
%% you need to change most of the paths in this block

%% preprocess for npx 2.0 probes : run kilosort, extract sync and lfp channels
function [channel_spikeCounts,channel_id] = bd_histologyRec(animal, date, site, shank, rerunHistCount)
animal = 'JF090';
date = '2022-11-17';
site = 1;
shank = 0;
rerunHistCount = 1;

%% convert raw data from .cbin to .bin and get channel map
myPaths;
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, [], 'ephys_histology', site, [], shank);
if size(ephysAPfile, 2) > 1 && iscell(ephysAPfile) %keep only ap
    ephysAPfile = ephysAPfile{1};
end

rootZ = fileparts(ephysAPfile);
metaFile = strrep(ephysAPfile, '.bin', '.meta');
[~, channelMapIMRO] = bc_readSpikeGLXMetaFile(metaFile);

%% get channel map file
chanMapFile = JF_imroToChannelMapLoc(channelMapIMRO);

%% kilosort
rootH = [extraHDPath, '/data_temp/'];
pathToYourConfigFile = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles'];
chanMapFilePath = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles', chanMapFile];
saveFile = [rootZ, filesep, '..', filesep, 'kilosort2', filesep, 'site', num2str(site), filesep];
nChannels = 385;

saveDir = dir(saveFile);
if isempty(saveDir)
    mkdir(saveFile)
end
%pyKS_saveFile = [toot, filesep, 'pykilosort' filesep, 'site', num2str(site), filesep, 'output', filesep];
histologyCountExists = dir(fullfile(saveFile, 'channels.spike_count.npy'));
if isempty(histologyCountExists) || rerunHistCount

    %addpath(genpath('D:\GitHub\KiloSort2')) % path to kilosort folder
    %addpath('D:\GitHub\npy-matlab') % for converting to Phy


    ops.trange = [0, Inf]; % time range to sort
    ops.NchanTOT = nChannels; % total number of channels in your recording

    run('/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/configKS2.m')
    ops.fproc = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
    ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);

    %% this block runs all the steps of the algorithm
    fprintf('Looking for data inside %s \n', rootZ)
    setenv('MW_NVCC_PATH', '/usr/lib/nvidia-cuda-toolkit/bin')
    % is there a channel map file in this folder?
    %fs = dir(fullfile(rootZ, 'chan*.mat'));
    % if ~isempty(fs)
    %     ops.chanMap = fullfile(rootZ, fs(1).name);
    % end

    % find the binary file
    fs = [dir(fullfile(rootZ, '*.bin')), dir(fullfile(rootZ, '*.dat'))];
    for iFolder = 1:size(fs, 1)
        useThis(iFolder) = contains(fs(iFolder).name, 'tcat');
    end
    if ~isempty(find(useThis))
        ops.fbinary = fullfile(rootZ, fs(useThis).name);

    else
        ops.fbinary = fullfile(rootZ, fs(1).name);

    end

    % preprocess data to create temp_wh.dat
    [igood, ich] = bd_preprocessDataSub(ops);
    [channel_spikeCounts,channel_id] = groupcounts(ich) ;
    channel_spikeCounts = gather(channel_spikeCounts); %gpu array to regular;
    channel_id = gather(channel_id); %gpu array to regular;

    writeNPY(channel_spikeCounts, fullfile(saveFile, 'channels.spike_count.npy'))
    writeNPY(channel_id, fullfile(saveFile, 'channels.id.npy'))
    %writeNPY()
else
    channel_spikeCounts = readNPY(fullfile(saveFile, 'channels.spike_count.npy'));
    channel_id = readNPY(fullfile(saveFile, 'channels.id.npy'));
end

end