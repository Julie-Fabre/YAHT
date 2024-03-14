function coordinates = ya_convert_allen_to_paxinos(values_toConvert, allenOrBrainreg)
% values to convert: matrix of ML x AP x DV values to convert 
% allenOrBrainreg: string, either 'allen' or 'brainreg'

% (scaling "Toronto MRI transform", reflect AP/ML, convert 25/10um to 1mm)
if strcmp(allenOrBrainreg, 'allen')
    scale = [0.952, -1.031, 0.885] ./ 100; % [ML, AP, DV]
elseif  strcmp(allenOrBrainreg, 'brainreg')
    scale = [0.952, -1.031, 0.885] ./ 250; % [ML, AP, DV]
else
    error('unrecognized value format')
end

% Scale the coordinates
AP_scaled = values_toConvert(2,:) * scale(2);
ML_scaled = values_toConvert(1,:) * scale(1);
DV_scaled = values_toConvert(3,:) * scale(3);

% Tilt angle
theta_degrees = 15; % Tilt by 15 degrees - this is my best approximation
theta_radians = deg2rad(theta_degrees);

% Rotation matrix for tilting around ML axis
R = [cos(theta_radians) 0 sin(theta_radians); 0 1 0; -sin(theta_radians) 0 cos(theta_radians)];

% Apply rotation to the scaled AP and DV (ignoring ML in rotation)
coordinates_scaled = [AP_scaled; ML_scaled; DV_scaled];
coordinates = R * coordinates_scaled;

end