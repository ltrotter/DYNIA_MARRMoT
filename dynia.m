function [cd_gradient,...
          cd_gradient_breaks,...
          info_content] = ...
                          dynia(model, n, Qobs, of_name, window, file_prefix)

if nargin < 6 || isempty(file_prefix); file_prefix = 'DYNIA'; end
file_theta = [file_prefix, '_theta_samples.csv'];
file_Qsim  = [file_prefix, '_Q_sim.csv'];
file_perf  = [file_prefix, '_OF_value_w', int2str(window), '.csv'];

% first check if any chunk was already run (i.e. if this is a rerun
% because the system ran out of time) and subtract from the number that
% still needs running
if isfile(file_theta)
    n_to_do = n - numel(textread(file_theta,'%1c%*[^\n]'));
else
    n_to_do = n;
end

chunks  = round(n_to_do/1000);    % divide it into chuncks of rougly 1000 points
n_chunk = round(n_to_do/chunks);  % number of points per chunck

% for each chunk
while chunks > 0
    % create Monte-Carlo sample of parameter sets
    theta_sample_chunk = unif_sample_par(model,n_chunk);
    
    % run the model with all of the samples created
    Qsim_chunk = run_with_par_sample(model, theta_sample_chunk);

    % calculate performance over a moving window of width window
    perf_over_time_chunk = calc_of_moving_window(Qsim_chunk, Qobs, window, of_name);

    % write both to file (appending to make sure you don't lose the values 
    % from the previous chunks), so that it can be retrieved afterwards
    writematrix(theta_sample_chunk', file_theta, "WriteMode", "append");
    writematrix(Qsim_chunk', file_Qsim, "WriteMode", "append");
    writematrix(perf_over_time_chunk', file_perf, "WriteMode", "append");

    chunks = chunks - 1;
end

% after all simulations are finished; open the three final files and
% continue
theta_sample = readmatrix(file_theta)';
Qsim = readmatrix(file_Qsim)';
perf_over_time = readmatrix(file_perf_over_time)';

% from this point on, timesteps with missing Qsim (where we couldn't
% calculate the objective function) can just be ignored, as long as we keep
% track of their position to re-add them in the end as missing rows
idx_non_missing = find(~isnan(perf_over_time(:,1)));
perf_over_time_non_missing = perf_over_time(idx_non_missing,:);

% keep top 10% of performances, which_top tells you the idx of the
% parameters used
[perf_top_10, which_top_10] = maxk(perf_over_time_non_missing, round(n/10),2);
% make sure all performances are >=0 by adding the minimum performance
perf_top_10_pos = perf_top_10 - min(min(perf_top_10,[],'all'),0);

% create empty container for best thetas for each timestep
theta_top_10 = zeros(numel(idx_non_missing), round(n/10), model.numParams);
% populate with the best sets of theta for each timestep
for t=1:size(theta_top_10,1)
    theta_top_10(t,:,:) = theta_sample(:,which_top_10(t,:))';
end

% calculate cumulative distributions
[cd_top_10, theta_top_10_sorted] = calc_cd_multidim(theta_top_10, perf_top_10_pos);

% calculate information content for each parameter at each timestep
info_content_non_missing = calc_info_content(cd_top_10, theta_top_10_sorted);

% calculate gradients of the (binned) cumulative distribution
[cd_gradient_breaks, cd_gradient_non_missing] = calc_cd_gradient(cd_top_10, theta_top_10_sorted);

%re-add in the timesteps with missing data
info_content = NaN(numel(Qobs), model.numParams);
info_content(idx_non_missing,:) = info_content_non_missing;

cd_gradient = NaN(numel(Qobs), size(cd_gradient_non_missing,2), model.numParams);
cd_gradient(idx_non_missing,:,:) = cd_gradient_non_missing;
