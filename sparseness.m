function s = sparseness(spike_counts)
% s = SPARSENESS(spike_counts)
% Calculates sparseness of a vector of positive integers 
% (baseline-normalized spike counts) following Vogels 1999.
% Close to 0: very sparse. 1: totally uniform.
if ~isvector(spike_counts)
    error('Spike counts should be a 1D vector')
end

n = length(spike_counts);
numerator = (sum(spike_counts/n))^2;
denominator = sum((spike_counts.^2)/n);
s = numerator / denominator;