function [atlasLocation, imgToRegister, imgToTransform, outputDir] = ...
    bd_getLocations(brainglobeLocation, brainsawPath, channelColToRegister, channelColToTransform,...
    atlasType, atlasSpecies, atlasResolution_um)


atlasLocation = dir([brainglobeLocation, atlasType, '_', ...
atlasSpecies, '_', num2str(atlasResolution_um), 'um*']); % atlas location

imgToRegister = dir([brainsawPath, filesep, '*', channelColToRegister, '*.tif*']);
imgToTransform = dir([brainsawPath, filesep, '*', channelColToTransform, '*.tif*']);

imgToTransform = imgToTransform(1);
imgToRegister = imgToRegister(1);
outputDir = [imgToRegister(1).folder, filesep, 'brainReg'];
end