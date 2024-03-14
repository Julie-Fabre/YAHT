function coordinates_bregma_mm = ya_convert_allen_to_paxinos(values_toConvert, allenOrBrainreg)
% values to convert: matrix of ML x AP x DV values to convert 
% allenOrBrainreg: string, either 'allen' or 'brainreg'

% (scaling "Toronto MRI transform", reflect AP/ML, convert 25/10um to 1mm)
if strcmp(allenOrBrainreg, 'allen')
    scale = [0.952, 1.031, 0.885] ./ 100; % [ML, AP, DV]
    bregma = [570, 540, 0]./100; % best estimate, [ML, AP, DV]
elseif  strcmp(allenOrBrainreg, 'brainreg')
    scale = [0.952, 1.031, 0.885] ./ 250; % [ML, AP, DV]
    bregma = [228, 216, 0]./ 250; % best estimate, [ML, AP, DV]
else
    error('unrecognized value format')
end

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

% Combine transformations: Note the order of multiplication
ccf_bregma_tform_matrix = ccf_translation_tform * ccf_scale_tform * ccf_rotation_tform;

% If using the transformation matrix directly in MATLAB, e.g., with affine3d
ccf_bregma_tform = affine3d(ccf_bregma_tform_matrix);

% Adjusting values relative to bregma
coordinates_bregma_mm = transformPointsForward(ccf_bregma_tform, values_toConvert);



end