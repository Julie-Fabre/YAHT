function bd_applyBrainReg(chanToTransform,outputDir)
% The deformation field files represent the mapping from one coordinate space (in brainreg’s case, 
% the downsampled, reoriented raw data space, i.e. the one in downsampled.tiff) to another (atlas space). 
% They are used as a lookup, and there are three files (0, 1, 2) to tell you the coordinate in atlas 
% space for the three spatial dimensions.
% 
% If you wanted to find out where a coordinate (say (100, 100, 100) in downsampled.tiff was in atlas space, you would:
% 
%     Find the first coordinate (this is based on the atlas, but typically anterior-posterior) by 
% finding the value of deformation_field_0.tiff at `(100,100,100).
%     Find the second coordinate (this is based on the atlas, but typically superior-inferior) by 
% finding the value of deformation_field_1.tiff at `(100,100,).100).
%     Find the third coordinate (this is based on the atlas, but typically right-left) by finding 
% the value of deformation_field_2.tiff at `(100,100,100).
% 
% You would also need to convert from mm (the scale of the deformation fields) to um (the scale of
% the atlas). This process is used in cellfinder here.
% 
% You could get a rough map of “transformation magnitude” by performing this process at every coordinate,
% and finding the geometric mean of the three values.




end