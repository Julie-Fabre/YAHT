
%% run brainreg
function [status, result] = bd_brainreg(channelToRegister, outputDir, orientationType)

%brainreg command . more details here: https://docs.brainglobe.info/brainreg/user-guide
CMD = sprintf('brainreg %s %s -v 25 25 25 --orientation %s %s ', ...
    channelToRegister, outputDir, orientationType);

%store a copy of the command to the directory
cmdFid = fopen(fullfile(outputDir, 'CMD'), 'w');
fprintf(cmdFid, '%s\n', CMD);
fclose(cmdFid);


% Run the command and report back if it failed
fprintf('Running: %s\n', CMD)
[status, result] = system(CMD);

end