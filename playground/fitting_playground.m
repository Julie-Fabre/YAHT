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

% Creating fittype objects for smooth spline
ft = fittype('smoothingspline');

% Fitting the curves
fitresult_x = fit(t, x, ft);
fitresult_y = fit(t, y, ft);
fitresult_z = fit(t, z, ft);

t_new = linspace(min(t), max(t), 1000)'; % Creating a new set of t values for smoother curve
x_new = fitresult_x(t_new); % Evaluating x at the new t values
y_new = fitresult_y(t_new); % Evaluating y at the new t values
z_new = fitresult_z(t_new); % Evaluating z at the new t values

% Plotting the original points and the fitted curve
figure;
scatter3(x, y, z, 'b'); % Original points in blue
hold on;
plot3(x_new, y_new, z_new, 'r'); % Fitted curve in red
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
hold off;
