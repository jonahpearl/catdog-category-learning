

% Generate some random data
data = randn(100, 1);

% Use MATLAB's built-in Gaussian kernel
[f, xi_builtin] = ksdensity(data);

% Use the custom Gaussian kernel
[fc, xi] = ksdensity(data, 'Kernel', @gaussian_kernel);

[fa, ~] = ksdensity(data, 'Kernel', @asymmetric_gaussian_kernel);

% Plotting the results to compare
figure;
plot(xi_builtin, f, 'r-', 'LineWidth', 2);
hold on;
plot(xi, fc, 'b--', 'LineWidth', 2);
plot(xi, fa, 'g--', 'LineWidth', 2);
legend('Built-in Gaussian', 'Custom Gaussian', 'Asymm Gaussian');
title('Comparison of Density Estimates');
xlabel('Data Points');
ylabel('Density');

%%
function y = gaussian_kernel(x)
    sigma = 1; % Standard deviation of the Gaussian kernel
    y = exp(-0.5 * (x / sigma).^2) / (sigma * sqrt(2 * pi));
end

function y = asymmetric_gaussian_kernel(x)
    sigma1 = 15; % SD for the causal side (positive x)
    sigma2 = 5;  % SD for the acausal side (negative x)
    
    % Allocate output array
    y = zeros(size(x));
    
    % Indices for positive and negative x
    positive = x >= 0;
    negative = x < 0;
    
    % Calculate kernel values
    y(positive) = exp(-0.5 * (x(positive) / sigma1).^2) / (sigma1 * sqrt(2 * pi));
    y(negative) = exp(-0.5 * (x(negative) / sigma2).^2) / (sigma2 * sqrt(2 * pi));
end
