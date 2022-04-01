function [cd_gradient,...
          cd_gradient_breaks,...
          info_content,...
          OF_idx] = ...
                          dynia(model, Qobs, options)

defaultopt = struct('repeats', 2e5, ...
                    'of_name', 'of_KGE', ...
                    'window_size', 31, ...
                    'window_step', 1,...
                    'file_prefix', 'DYNIA', ...
                    'chunk_size', 1000,...
                    'precision_Q', 4, ...
                    'precision_OF', 4);
defaultopt.of_args = cell(0);

if nargin < 3 || isempty(options); options = struct(); end

% get options
n = optimget(options, 'repeats', defaultopt, 'fast');
of_name = optimget(options, 'of_name', defaultopt, 'fast');
window = optimget(options, 'window_size', defaultopt, 'fast');
step = optimget(options, 'window_step', defaultopt, 'fast');
file_prefix = optimget(options, 'file_prefix', defaultopt, 'fast');
c_size = optimget(options, 'chunk_size', defaultopt, 'fast');
precision_Q = optimget(options, 'precision_Q', defaultopt, 'fast');
precision_OF = optimget(options, 'precision_OF', defaultopt, 'fast');
of_args = optimget(options, 'of_args', defaultopt, 'fast');

file_theta = [file_prefix, '_theta_samples.csv'];
file_Qsim  = [file_prefix, '_Q_sim'];
file_perf  = [file_prefix, '_OF_value'];
file_log   = [file_prefix, '.mat'];

% calculate the total number of chunks 

% first check if any chunk was already run (i.e. if this is a rerun
% because the system ran out of time) and make sure the options used were
% the same, then subtract the ones already ran to the total number of
% repeats.
% if this is the first time you run this, write the options used at the top
% of the file, so that they can be checked afterwards

if any(~isfile(file_log))
    n_done  = 0;
    OF_idx = (floor(window/2)+1):step:(numel(Qobs)-floor(window/2)-1);
    save(file_log, "of_name", "window", "step", "of_args", "OF_idx", "n_done");
else
    disp([file_log ' found: some options will be loaded.'])
    load(file_log, "of_name", "window", "step", "of_args", "OF_idx", "n_done")
end

n_to_do = max(0,n-n_done);
chunks  = round(n_to_do/c_size);  % divide it into chuncks of rougly c_size points
n_chunk = ceil(n_to_do/chunks);  % actual number of points per chunck

% for each chunk
while chunks > 0
    % create Monte-Carlo sample of parameter sets
    theta_sample_chunk = unif_sample_par(model,n_chunk);

    % run the model with all of the samples created
    Qsim_chunk = run_with_par_sample(model, theta_sample_chunk);

    % calculate performance over a moving window of width window
    perf_over_time_chunk = calc_of_moving_window(Qsim_chunk, Qobs, window, step, of_name, of_args{:});

    % write both to file (appending to make sure you don't lose the values
    % from the previous chunks), so that it can be retrieved afterwards
    writematrix(theta_sample_chunk', file_theta, "WriteMode", "append");
    writematrix(round(Qsim_chunk, precision_Q)', [file_Qsim,'_',num2str(chunks,'%03d' ),'.csv'], "WriteMode", "append");
    writematrix(round(perf_over_time_chunk, precision_OF)', [file_perf,'_',num2str(chunks,'%03d'),'.csv'], "WriteMode", "append");

    chunks = chunks - 1;
    n_done = n_done + n_chunk;
    save(file_log, "n_done", "-append");
end

% after all simulations are finished; open the three final files and
% continue
theta_sample = readmatrix(file_theta)';
perf_ds = datastore([file_perf, '_*.csv']);
perf_over_time = table2array(perf_ds.readall())';
%perf_over_time = readmatrix(file_perf, 'NumHeaderLines',2)';

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
info_content = NaN(numel(OF_idx), model.numParams);
info_content(idx_non_missing,:) = info_content_non_missing;

cd_gradient = NaN(numel(OF_idx), size(cd_gradient_non_missing,2), model.numParams);
cd_gradient(idx_non_missing,:,:) = cd_gradient_non_missing;

save(file_log, "cd_gradient", "cd_gradient_breaks", "info_content", "-append");