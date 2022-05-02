function [cd, x_sorted] = calc_cd_multidim(x,y)
% x is expected to be a vector of n parameter sets at each timestep
% its dimensions are T,n,P where P is the number of parameters and T the
% number of timesteps

% y is expected to be a vector of the n best performance values at each 
% timesteps (the performance associated with each parameter set in x)
% its dimensions are T,n

% sort each parameter at each timestep from smallest to largest
[x_sorted, idx_x_sorted] = sort(x, 2);

y_normal = y ./ sum(y,2, 'omitnan');

% create empty container for the cum distributions
cd = zeros(size(x));
% for each timestep
for t=1:size(x, 1)
    % for each parameter
    for p=1:size(x,3)
        % calc CD and populate the cd_top_10 vector 
        cd(t,:,p) = cumsum(y_normal(t,idx_x_sorted(t,:,p)));
    end
end