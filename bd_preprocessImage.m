function imgToRegister_processed = bd_preprocessImage(imagePath, outputDir)

image_toProcess = loadtiff(imagePath);

%% remove outlier values 
lowpass_value = min(prctile(image_toProcess(:), 5), 10);
highpass_value = max(prctile(image_toProcess(:), 95), 245);

image_toProcess(image_toProcess < lowpass_value) = NaN;
image_toProcess(image_toProcess > highpass_value) = NaN;

%% despeckle 
window_size = [3, 3, 3];
filtered_image = medfilt3(image_toProcess, window_size);

saveastiff(filtered_image, [imagePath(1:end-4), '_filtered.tiff'])
imgToRegister_processed = dir([imagePath(1:end-4), '_filtered.tiff']);
% figure(); imagesc(filtered_image(:,:,250)); colormap("gray")
% figure(); imagesc(image_toProcess(:,:,250)); colormap("gray")
% option 2 
% K = wiener2(J,[5 5]);
end