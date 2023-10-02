% t = (1:5)';
% x = [1; 1.8; 2.4; 3.5; 4.2];
% y = [2; 2.5; 3.2; 3.9; 4.8];
% z = [3; 3.3; 3.8; 4.4; 5.1];
% 
% % Fit using pchip
% pp_x = pchip(t, x);
% pp_y = pchip(t, y);
% pp_z = pchip(t, z);
% 
% % Compute the gradients (tangents) at the last data point
% pp_deriv_x = fnder(pp_x);
% pp_deriv_y = fnder(pp_y);
% pp_deriv_z = fnder(pp_z);
% 
% grad_x = ppval(pp_deriv_x, max(t));
% grad_y = ppval(pp_deriv_y, max(t));
% grad_z = ppval(pp_deriv_z, max(t));
% 
% % Determine the range for extrapolation and steps
% t_extended = (6:8)';
% delta_t = t_extended(2) - t_extended(1);  % Step difference
% 
% % Extrapolate using the gradients
% x_extended = x(end) + cumsum(grad_x * delta_t);
% y_extended = y(end) + cumsum(grad_y * delta_t);
% z_extended = z(end) + cumsum(grad_z * delta_t);
% 
% % Append the data
% t_total = [t; t_extended];
% x_total = [x; x_extended];
% y_total = [y; y_extended];
% z_total = [z; z_extended];
% 
% % Plot
% figure;
% plot3(x, y, z, 'bo-', x_extended, y_extended, z_extended, 'ro-');
% legend('Original data', 'Extrapolated data');
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% 
% 
% This method ensures that the extrapolated data points continue in the direction given by the tangent vector at the last data point. Adjusting the delta_t can give more or less aggressive extrapolations.
% 
% 
% % Sample data
% t = (1:5)';
% x = [1; 1.8; 2.9; 3.5; 4.1];
% 
% % Smooth the data using csaps (with 0.9 as the smoothing parameter)
% p = 0.9;
% sp_x = csaps(t, x, p);
% 
% % Convert the smoothed spline to points for enforcing monotonicity
% t_dense = linspace(min(t), max(t), 1000);
% x_smooth = fnval(sp_x, t_dense);
% 
% % Enforce monotonicity on the smoothed curve
% for i = 2:length(x_smooth)
%     x_smooth(i) = max(x_smooth(i), x_smooth(i-1));  % For increasing monotonicity
%     % x_smooth(i) = min(x_smooth(i), x_smooth(i-1)); % Uncomment for decreasing monotonicity
% end
% 
% % Extrapolate
% t_extended = linspace(max(t)+0.1, max(t)+2, 100); % for example, extending 2 units
% x_extended = x_smooth(end) + cumsum((x_smooth(end)-x_smooth(end-1)) * diff([t_dense(end), t_extended]));
% 
% % Combine data
% t_total = [t_dense, t_extended];
% x_total = [x_smooth, x_extended];
% 
% % Plot
% figure;
% plot(t, x, 'bo', t_total, x_total, 'r-');
% legend('Original Data', 'Smoothed Monotonic Curve');
% grid on;
% 
% % Sample data
% t = (1:5)';
% x = [1; 1.8; 2.9; 3.5; 4.1];
% 
% % Smooth the data using csaps
% p = 0.9;
% sp_x = csaps(t, x, p);
% 
% % Convert the smoothed spline to points for enforcing monotonicity
% t_dense = linspace(min(t), max(t), 1000);
% x_smooth = fnval(sp_x, t_dense);
% 
% % Enforce monotonicity on the smoothed curve
% for i = 2:length(x_smooth)
%     x_smooth(i) = max(x_smooth(i), x_smooth(i-1));  % For increasing monotonicity
% end
% 
% % Improved Extrapolation using an average gradient
% segments_to_consider = 50;  % last 50 segments, for example
% average_gradient = mean(diff(x_smooth(end-segments_to_consider:end)));
% 
% t_extended = linspace(max(t)+0.1, max(t)+2, 100);
% x_extended = x_smooth(end) + cumsum(average_gradient * diff([t_dense(end), t_extended]));
% 
% % Combine data
% t_total = [t_dense, t_extended];
% x_total = [x_smooth, x_extended];
% 
% % Plot
% figure;
% plot(t, x, 'bo', t_total, x_total, 'r-');
% legend('Original Data', 'Smoothed Monotonic Curve');
% grid on;
% 
% 
% % Fit Gaussian Process Regression
% gprMdl = fitrgp(t,x);
% 
% % Predict, including extrapolation
% t_extended = linspace(min(t), max(t)+2, 1000)'; % extending 2 units for example
% x_predicted = predict(gprMdl, t_extended);
% 
% % Plot
% figure;
% plot(t, x, 'bo', t_extended, x_predicted, 'r-');
% legend('Original Data', 'GPR Fit');
% grid on;
% 
% 
% % Sample data
% t = (1:5)';
% x = [1; 1.8; 2.9; 3.5; 4.1];
% 
% % Piecewise Linear Interpolation
% t_dense = linspace(min(t), max(t)+2, 1000)';
% x_interpolated = interp1(t, x, t_dense, 'linear');
% 
% % Apply a Smoothing Filter
% windowSize = 50;
% x_smoothed = movmean(x_interpolated, windowSize);
% 
% % Enforce Monotonicity
% for i = 2:length(x_smoothed)
%     x_smoothed(i) = max(x_smoothed(i), x_smoothed(i-1)); % For increasing monotonicity
% end
% 
% % Plot
% figure;
% plot(t, x, 'bo', t_dense, x_smoothed, 'r-');
% legend('Original Data', 'Smoothed Monotonic Curve');
% grid on;
% 
% % Sample data
% t = (1:5)';
% x = [1; 1.8; 2.9; 3.5; 4.1];
% 
% % Calculate differences
% dx = diff(x);
% dt = diff(t);
% 
% % Fit a smoothing spline to the differences
% ft = fittype('smoothingspline');
% opts = fitoptions('Method', 'SmoothingSpline');
% opts.SmoothingParam = 0.01;  % Adjust this value for different smoothing levels
% fitresult = fit(t(1:end-1), dx, ft, opts);
% 
% % Integrate the fitted differences to get the curve
% x_fit = [x(1); x(1) + cumsum(feval(fitresult, t(1:end-1)))];
% 
% % Define the extended time range for extrapolation
% t_extended = (max(t)+1):(max(t)+5);  % For example, extending 5 units
% dt_extended = diff(t_extended);
% 
% % Predict the differences for the extended range using feval with fitresult
% dx_extended = feval(fitresult, t_extended(1:end)');
% 
% % Combine the original and extrapolated differences
% total_dx = [dx; dx_extended];
% 
% % Integrate the combined differences to get the total curve
% x_total = [x(1); x(1) + cumsum(total_dx)];
% 
% % Combine original and extended time values
% t_total = [t; t_extended'];
% 
% % Plot
% figure;
% plot(t, x, 'bo', t_total, x_total, 'r-');
% legend('Original Data', 'Monotonic Smoothed Curve with Extrapolation');
% grid on;
% 
% 
% % Sample data
% t = (1:5)';
% x = [1; 1.8; 2.9; 3.5; 4.1];
% 
% % Calculate differences
% dx = diff(x);
% 
% % Fit a smoothing spline to the differences
% ft = fittype('smoothingspline');
% opts = fitoptions('Method', 'SmoothingSpline');
% opts.SmoothingParam = 0.01;  % Adjust this value for different smoothing levels
% fitresult = fit(t(1:end-1), dx, ft, opts);
% 
% % Integrate the fitted differences to get the curve
% x_fit = [x(1); x(1) + cumsum(feval(fitresult, t(1:end-1)))];
% 
% % Define the extended time range for extrapolation
% t_extended = (max(t)+1):(max(t)+5);  % For example, extending 5 units
% 
% % Use the gradient of the latter segment of the smoothed data as trend direction
% trend_direction = mean(dx(end-2:end));  % Average of the last 3 differences
% dx_extended = repmat(trend_direction, length(t_extended), 1);
% 
% % Combine the original and extrapolated differences
% total_dx = [dx; dx_extended];
% 
% % Integrate the combined differences to get the total curve
% x_total = [x(1); x(1) + cumsum(total_dx)];
% 
% % Combine original and extended time values
% t_total = [t; t_extended'];
% 
% % Plot
% figure;
% plot(t, x, 'bo', t_total, x_total, 'r-');
% legend('Original Data', 'Monotonic Smoothed Curve with Extrapolation');
% grid on;


% Sample data
t = (1:5)';
x = [1; 1.8; 2.9; 3.5; 4.1];
y = [2; 2.5; 3.2; 3.7; 4.3];
z = [5; 5.2; 5.5; 6.0; 6.8];



x_total = process_3d_coord(t, x);
y_total = process_3d_coord(t, y);
z_total = process_3d_coord(t, z);



% Plot
figure;
plot3(x, y, z, 'bo', x_total, y_total, z_total, 'r-');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('Original Data', 'Monotonic Smoothed Curve with Extrapolation');
grid on;


% Sample 3D data
t = (1:5)';
x = [1; 1.8; 2.9; 3.5; 4.1];
y = [2; 2.5; 3.2; 3.7; 4.3];
z = [5; 5.2; 5.5; 6.0; 6.8];



% Obtain a smooth spline fit and then enforce monotonicity
smoothing_parameter = 0.9; % You can adjust this
t_extended = linspace(min(t), max(t)+5, 100);
x_smooth = csaps(t, x, smoothing_parameter, t_extended);
x_monotonic = enforce_monotonic(x_smooth);
y_smooth = csaps(t, y, smoothing_parameter, t_extended);
y_monotonic = enforce_monotonic(y_smooth);
z_smooth = csaps(t, z, smoothing_parameter, t_extended);
z_monotonic = enforce_monotonic(z_smooth);

% Plot
figure;
plot3(x, y, z, 'bo', x_monotonic, y_monotonic, z_monotonic, 'r-');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('Original Data', 'Monotonic Smoothed Curve');
grid on;

% Custom function to enforce monotonicity after smoothing

function out = enforce_monotonic(in)
    for i = 2:length(in)
        if in(i) < in(i-1)  % or use > for enforcing increasing monotonicity
            in(i) = in(i-1);
        end
    end
    out = in;
end

% Function to handle the monotonic smoothed curve with extrapolation
function [out_coord] = process_coord(t, coord)
    dx = diff(coord);
    ft = fittype('smoothingspline');
    opts = fitoptions('Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.01;
    fitresult = fit(t(1:end-1), dx, ft, opts);
    coord_fit = [coord(1); coord(1) + cumsum(feval(fitresult, t(1:end-1)))];
    t_extended = (max(t)+1):(max(t)+5);
    trend_direction = mean(dx(end-2:end));
    dx_extended = repmat(trend_direction, length(t_extended), 1);
    total_dx = [dx; dx_extended];
    out_coord = [coord(1); coord(1) + cumsum(total_dx)];
end

function out_coord = process_3d_coord(t, coord)
    % Fit a smoothing spline
    [fitresult, ~] = fit(t, coord, 'smoothingspline', 'SmoothingParam', 0.01);
    
    % Ensure monotonicity by taking the gradient direction of the end of the fitted data
    trend_direction = sign(mean(gradient(feval(fitresult, t(end-2:end)))));
    
    % Extend t for extrapolation
    t_extended = (max(t)+1):(max(t)+5);
    coord_extended = feval(fitresult, t_extended);
    
    % Enforce monotonicity on extrapolated part
    coord_extended = coord(end) + trend_direction * cumsum(abs(diff([coord(end); coord_extended])));
    
    out_coord = [feval(fitresult, t); coord_extended(2:end)];
end