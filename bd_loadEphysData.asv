function [spike_times, spike_templates, template_depths, spike_depths, ...
    template_xdepths, spike_xdepths, channelPositions, goodChannels] = bd_loadEphysData(ephys_path, ephys_sample_rate)
% JF, Load ephys data (1-indexed)
% ------
% Inputs
% ------
% ephys_path: character array defining the path to your kilosorted output files

% ------
% Outputs
% ------


spike_templates_0idx = readNPY([ephys_path, filesep, 'spike_templates.npy']);
spike_templates = spike_templates_0idx + 1;
if exist(fullfile(ephys_path, 'spike_times_corrected.npy')) % When running pyKS stitched you need the 'aligned / corrected' spike times
    spikeTimes_samples = double(readNPY([ephys_path, filesep, 'spike_times_corrected.npy']));
    spikeTimes_datasets = double(readNPY([ephys_path, filesep, 'spike_datasets.npy'])) + 1; %  which dataset? (zero-indexed so +1)
else
    spikeTimes_samples = double(readNPY([ephys_path, filesep, 'spike_times.npy']));
    spikeTimes_datasets = ones(size(spikeTimes_samples));
end
spike_times = spikeTimes_samples ./ param.ephys_sample_rate;

channelPositions = readNPY([ephys_path, filesep, 'channel_positions.npy']);
goodChannels = readNPY([ephys_path, filesep, 'channel_map.npy']) + 1;

% Load and unwhiten templates
templateWaveforms_whitened = readNPY([ephys_path filesep 'templates.npy']);
winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
templates = zeros(size(templateWaveforms_whitened));
for t = 1:size(templates,1)
    templates(t,:,:) = squeeze(templateWaveforms_whitened(t,:,:))*winv;
end

% Get the waveform of all templates (channel with largest amplitude)
[~, max_site] = max(max(abs(templates), [], 2), [], 3);
templates_max = nan(size(templates, 1), size(templates, 2));
for curr_template = 1:size(templates, 1)
    templates_max(curr_template, :) = ...
        templates(curr_template, :, max_site(curr_template));
end

% Get depth of each template
% (get min-max range for each channel)
template_chan_amp = squeeze(range(templates, 2));
% (zero-out low amplitude channels)
template_chan_amp_thresh = max(template_chan_amp, [], 2) * 0.5;
template_chan_amp_overthresh = template_chan_amp .* (template_chan_amp >= template_chan_amp_thresh);
% (get center-of-mass on thresholded channel amplitudes)
template_depths = sum(template_chan_amp_overthresh.*channel_positions(:, 2)', 2) ./ sum(template_chan_amp_overthresh, 2);
template_xdepths = channel_positions(max_site, 1);

% Get the depth of each spike (templates are zero-indexed)
spike_depths = template_depths(spike_templates_0idx+1);
spike_xdepths = template_xdepths(spike_templates_0idx+1);


end