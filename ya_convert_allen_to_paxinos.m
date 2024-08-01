function coordinates_bregma_mm = ya_convert_allen_to_paxinos(values_toConvert, allenOrBrainreg)
% based on AP_histology
% values to convert: matrix of ML x AP x DV values to convert 
% allenOrBrainreg: string, either 'allen' or 'brainreg'

values_toConvert_ori = values_toConvert;
if strcmp(allenOrBrainreg, 'allen')
    scaling_factor = 100; % um to mm
    bregma = [570.5,520,44];% ML x AP x DV
elseif  strcmp(allenOrBrainreg, 'brainreg')
    scaling_factor = 100; % um to mm
    values_toConvert = values_toConvert_ori .* [2.5, 2.5, 2.5]; % 25um/slice to 10um/slice (AP)
    bregma = [570.5,520,44];% ML x AP x DV
else
    error('unrecognized value format')
end

scale = [0.952, -1.031, 0.885] ./ scaling_factor; % ML x AP x DV;  (scaling "Toronto MRI transform", reflect AP/ML, convert 25/10um to 1mm)

% Translation transformation matrix
ccf_translation_tform = eye(4) + [zeros(3,4); -bregma, 0];

% Scaling transformation matrix
ccf_scale_tform = eye(4) .* [scale, 1]';

% Rotation (tilt) adjustment
ap_rotation = 15; % Tilting the CCF 15 degrees nose-up - this is my best approximation

% Rotation transformation matrix
ccf_rotation_tform = ...
    [1 0 0 0; ...
     0 cosd(ap_rotation) -sind(ap_rotation) 0; ...
     0 sind(ap_rotation) cosd(ap_rotation) 0; ...
     0 0 0 1];

% Combine transformations: is this order of multiplication correct?
ccf_bregma_tform_matrix =  ccf_translation_tform * ccf_scale_tform * ccf_rotation_tform; 

% Full transf. matrix
ccf_bregma_tform = affine3d(ccf_bregma_tform_matrix);

% Adjusting values relative to bregma
coordinates_bregma_mm = transformPointsForward(ccf_bregma_tform, values_toConvert);



end