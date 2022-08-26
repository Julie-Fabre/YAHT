function bd_applyBrainReg(imgToTransform, atlasResolution_um, atlasSize, outputDir)
% QQ this doesn't work very well - need to trouble shoot more 

% from https://forum.image.sc/t/whats-in-the-deformation-files/66409:
% The deformation field files represent the mapping from one coordinate space (in brainregfs case, 
% the downsampled, reoriented raw data space, i.e. the one in downsampled.tiff) to another (atlas space). 
% They are used as a lookup, and there are three files (0, 1, 2) to tell you the coordinate in atlas 
% space for the three spatial dimensions.
% 
% If you wanted to find out where a coordinate (say (100, 100, 100) in downsampled.tiff was in atlas space, you would:
% 
%     Find the first coordinate (this is based on the atlas, but typically anterior-posterior) by 
% finding the value of deformation_field_0.tiff at `(100,100,100).
%     Find the second coordinate (this is based on the atlas, but typically superior-inferior) by 
% finding the value of deformation_field_1.tiff at `(100,100,100).
%     Find the third coordinate (this is based on the atlas, but typically right-left) by finding 
% the value of deformation_field_2.tiff at `(100,100,100).
% 
% You would also need to convert from mm (the scale of the deformation fields) to um (the scale of
% the atlas). This process is used in cellfinder here:
%             
%             def transform_points_downsampled_to_atlas_space(
%                 downsampled_points, atlas, deformation_field_paths, output_filename=None
%             ):
%                 field_scales = [int(1000 / resolution) for resolution in atlas.resolution]
%                 points = [[], [], []]
%                 for axis, deformation_field_path in enumerate(deformation_field_paths):
%                     deformation_field = tifffile.imread(deformation_field_path)
%                     for point in downsampled_points:
%                         point = [int(round(p)) for p in point]
%                         points[axis].append(
%                             int(
%                                 round(
%                                     field_scales[axis]
%                                     * deformation_field[point[0], point[1], point[2]]
%                                 )
%                             )
%                         )
%             % 
% You could get a rough map of “transformation magnitude” by performing this process at every coordinate,
% and finding the geometric mean of the three values.

imgTr = loadtiff(imgToTransform);
imgOr = loadtiff(['/home/netshare/znas-brainsaw/JF070_JF078/JF070/downsampled_stacks/025_micron/brainReg/downsampled_standard.tiff']);


deformation1 = loadtiff([outputDir, filesep, 'deformation_field_0.tiff']);
deformation2 = loadtiff([outputDir, filesep, 'deformation_field_1.tiff']);
deformation3 = loadtiff([outputDir, filesep, 'deformation_field_2.tiff']);

transformedImage = nan([atlasSize(3), atlasSize(1), atlasSize(2)]); 
%transformedImage = nan
for iPixelX = 1:size(imgTr,1) % QQ make this faster ...
    for iPixelY = 1:size(imgTr,2)
        for iPixelZ = 1:size(imgTr,3)
            locationV = round(1000/ atlasResolution_um * [deformation1(iPixelX, iPixelY, iPixelZ),...
    deformation2(iPixelX, iPixelY, iPixelZ), deformation3(iPixelX, iPixelY, iPixelZ)]);
            if ~any(locationV < 1) && ~any(locationV > [atlasSize(3), atlasSize(1), atlasSize(2)])
                transformedImage(locationV(1), locationV(2), locationV(3)) = imgTr(iPixelX, iPixelY, iPixelZ);
            end
        end
    end
end
saveFile = transformedImage;
%saveastiff(saveFile , [outputDir, filesep, 'downsampled_transformed_312p.tiff']); 

figure(); 
subplot(311); imagesc(squeeze(imgTr(:, :, 100))); axis equal;

subplot(312); imagesc(squeeze(transformedImage(100, :, :))); axis equal;
%subplot(412); imagesc(squeeze(transformedImage(:, 100, :))); axis equal;
%subplot(413); imagesc(squeeze(transformedImage(:, :, 100))); axis equal;
subplot(313); imagesc(squeeze(imgOr(:, :, 528-100))); axis equal;

end