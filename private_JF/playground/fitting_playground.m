%% fit probes playground 

load('/home/netshare/zaru/JF097/Histology/downsampled_stacks/025_micron/brainReg/probe_points.mat')


iProbe = 1;
slice_points = find(~cellfun(@isempty, probe_points(:,1)));

points_per_slice = arrayfun(@(x) size(probe_points{slice_points(x), iProbe},1), 1:length(slice_points));
n_points = sum(points_per_slice);

probe_points_full = zeros(n_points, 3);
thisCount = 1;
for iSlice = 1:size(slice_points,1)
    for iPoint = 1:points_per_slice(iSlice)
        probe_points_full(thisCount,:) = [slice_points(iSlice), probe_points{slice_points(iSlice), iProbe}(iPoint,:)];
        thisCount = thisCount + 1;
    end
end
x = probe_points_full(:,1);
y = probe_points_full(:,2);
z = probe_points_full(:,3);

t = (1:length(x))'; % Creating a parameter t based on the number of data points

% Creating fittype objects for smooth spline with smoothing parameter
ft = fittype('smoothingspline');

% Define a smoothing parameter
smoothing_param = 0.1; % Adjust this value between 0 and 1 to control the smoothness

% Fitting the curves with the smoothing parameter
fitresult_x = fit(t, x, ft, 'SmoothingParam', smoothing_param);
fitresult_y = fit(t, y, ft, 'SmoothingParam', smoothing_param);
fitresult_z = fit(t, z, ft, 'SmoothingParam', smoothing_param);

delta = 5;
% Step 3: Extending the Curve and Visualizing
% Extend the range of t for extrapolation
t_extended = [linspace(min(t)-delta, max(t)+delta, 1000)]'; % Adjust delta to control the extension amount

% Evaluate the extended spline at the new t values
x_new = fitresult_x(t_extended);
y_new = fitresult_y(t_extended);
z_new = fitresult_z(t_extended);

% Plot the original points and the extended spline
figure;
scatter3(x, y, z, 'b'); % Original points in blue
hold on;
plot3(x_new, y_new, z_new, 'r'); % Extended spline in red
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
hold off;


gui_data.spline_smoothness_factor
