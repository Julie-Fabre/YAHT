function [atlasLocation, imgToRegister, imgToTransform, outputDir] =...
    bd_getLocations(brainglobeLocation, brainsawPath, animal, channelColToRegister, channelColToTransform, atlasType, atlasSpecies, atlasResolution_um)
atlasLocation = dir([brainglobeLocation, atlasType, '_',...
    atlasSpecies, '_', num2str(atlasResolution_um), 'um*']); % atlas location
imgToRegister = dir([brainsawPath, '/*/', animal, '/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
imgToTransform = dir([brainsawPath, '/*/', animal, '/downsampled_stacks/025_micron/*', channelColToTransform, '*.tif*']);
if isempty(imgToRegister) || contains(imgToRegister.folder, 'Recycle')
    imgToRegister = dir(['/home/netshare/zaru/', animal, '/*istology/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
    imgToTransform = dir(['/home/netshare/zaru/', animal, '/*istology/downsampled_stacks/025_micron/*', channelColToTransform, '*.tif*']);
end
outputDir = [imgToRegister.folder, filesep, 'brainReg'];
end