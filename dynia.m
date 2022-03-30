function [cd_gradient, info_content, theta_sample, Qsim, perf_over_time] = ...
                                     dynia(model, n, Qobs, of_name, window)

% create Monte-Carlo sample of parameter sets
theta_sample = unif_sample_par(model,n);

% run the model with all of the samples created
Qsim = run_with_par_sample(model, theta_sample);

% calculate performance over a moving window of width window
perf_over_time = calc_of_moving_window(Qsim, Qobs, window, of_name);

% keep top 10% of performances, which_top tells you the idx of the
% parameters used
[perf_top_10, which_top_10] = maxk(perf_over_time, round(n/10),2);
% make sure all performances are >=0 by adding the minimum performance
perf_top_10_pos = perf_top_10 - min(min(perf_top_10,[],'all'),0);

% create empty container for best thetas for each timestep
theta_top_10 = zeros(numel(Qobs), round(n/10), model.numParams);
% populate with the best sets of theta for each timestep
for t=1:size(theta_top_10,1)
    theta_top_10(t,:,:) = theta_sample(:,which_top_10(t,:))';
end

% calculate cumulative distributions
[cd_top_10, theta_top_10_sorted] = calc_cd_multidim(theta_top_10, perf_top_10_pos);

% calculate information content for each parameter at each timestep
info_content = calc_info_content(cd_top_10, theta_top_10_sorted);

% calculate gradients of the (binned) cumulative distribution
cd_gradient = calc_cd_gradient(cd_top_10, theta_top_10_sorted);