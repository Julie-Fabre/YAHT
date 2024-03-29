
%% run brainreg
function [status, result] = ya_brainreg(channelToRegister, outputDir, probeTrackChannel, orientationType, atlas)

if atlas == 10 
    atlasString = 'allen_mouse_10um';
elseif atlas == 25 
    atlasString = 'allen_mouse_25um';
end

%brainreg command . more details here: https://docs.brainglobe.info/brainreg/user-guide
CMD = sprintf('brainreg %s %s --downsample %s -v 25 25 25 --orientation %s --atlas %s', ...
    channelToRegister, outputDir, probeTrackChannel, orientationType, atlasString);

%store a copy of the command to the directory
mkdir(outputDir)
cmdFid = fopen(fullfile(outputDir, 'CMD'), 'w');
fprintf(cmdFid, '%s\n', CMD);
fclose(cmdFid);


% Run the command and report back if it failed
fprintf('Running: %s\n', CMD)
[status, result] = system(CMD);

if status ~= 0 
    error('error: brain registration failed ')
end
% ouput file information :
% https://docs.brainglobe.info/brainreg/user-guide/output-files 
end