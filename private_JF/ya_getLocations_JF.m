function [atlasLocation, imgToRegister, imgToTransform, outputDir] =...
    ya_getLocations(brainglobeLocation, brainsawPath, animal, channelColToRegister, channelColToTransform, atlasType, atlasSpecies, atlasResolution_um)
atlasLocation = dir([brainglobeLocation, atlasType, '_',...
    atlasSpecies, '_', num2str(atlasResolution_um), 'um*']); % atlas location
% default in subjects 
imgToRegister = dir(['/home/netshare/zaru/', animal, '/*istology/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
imgToTransform = dir(['/home/netshare/zaru/', animal, '/*istology/downsampled_stacks/025_micron/*', channelColToTransform, '*.tif*']);


if isempty(imgToRegister) 
    imgToRegister = dir(['/home/netshare/zinu/', animal, '/*istology/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
    imgToTransform = dir(['/home/netshare/zinu/', animal, '/*istology/downsampled_stacks/025_micron/*', channelColToTransform, '*.tif*']);
end
if isempty(imgToRegister) 
    imgToRegister = dir(['/home/netshare/znas/', animal, '/*istology/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
    imgToTransform = dir(['/home/netshare/znas/', animal, '/*istology/downsampled_stacks/025_micron/*', channelColToTransform, '*.tif*']);
end
% otherwise in znas-brainsaw-display warning 
if isempty(imgToRegister) || contains(imgToRegister(1).folder, 'Recycle')
    imgToRegister = dir([brainsawPath, '/*/', animal, '/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
    imgToTransform = dir([brainsawPath, '/*/', animal, '/downsampled_stacks/025_micron/*', channelColToTransform, '*.tif*']);
end
if isempty(imgToRegister)
    imgToRegister = dir([brainsawPath, animal, '/downsampled_stacks/025_micron/*', channelColToRegister, '*.tif*']);
    imgToTransform = dir([brainsawPath, animal, '/downsampled_stacks/025_micron/*', channelColToTransform, '*.tif*']);
end
imgToTransform = imgToTransform(1);
imgToRegister = imgToRegister(1);
outputDir = [imgToRegister(1).folder, filesep, 'brainReg'];
end