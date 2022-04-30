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
    n_done = 0;
    last_fid = 0;
    OF_idx = (floor(window/2)+1):step:(numel(Qobs)-floor(window/2)-1);     % this is a raw estimate, after the calculation of the OF, it will be updated
    save(file_log, "of_name", "window", "step", "of_args", "OF_idx", "n_done", "last_fid");
else
    disp([file_log ' found: some options will be loaded.'])
    load(file_log, "of_name", "window", "step", "of_args", "OF_idx", "n_done", "last_fid");
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
    [OF_idx, perf_over_time_chunk] = calc_of_moving_window(Qsim_chunk, Qobs, window, step, of_name, precision_Q+1, of_args{:});

    % write both to file (appending to make sure you don't lose the values
    % from the previous chunks), so that it can be retrieved afterwards
    writematrix(theta_sample_chunk', file_theta, "WriteMode", "append");
    last_fid = last_fid + 1;
    writematrix(round(Qsim_chunk, precision_Q)', [file_Qsim,'_',num2str(last_fid,'%03d' ),'.csv']);
    writematrix(round(perf_over_time_chunk, precision_OF)', [file_perf,'_',num2str(last_fid,'%03d'),'.csv']);

    chunks = chunks - 1;
    n_done = n_done + n_chunk;
    save(file_log, "n_done", "last_fid", "-append");
    if chunks == 1; save(file_log, "OF_idx", "-append"); end
end

% after all simulations are finished; open the final files and
% continue
theta_sample = readmatrix(file_theta)';
perf_ds = datastore([file_perf, '_*.csv']);

% create empty containers for top10 performance and parameters at every
% timestep
which_top_10 = zeros(numel(OF_idx), round(n/10));
perf_top_10 = zeros(numel(OF_idx), round(n/10));

% divide this in chunks too otherwise it is too large:
e = 0;
c_size2 = round(c_size/(n/numel(Qobs))); 
while e < numel(OF_idx)
    disp(["extracting performance: ", int2str(s) ":" int2str(e), "/ of", int2str(numel(OF_idx))])
    s = e + 1; e = min(s + c_size2 -1, numel(OF_idx));
    % select variable for each chunk
    perf_ds.SelectedVariableNames = arrayfun(@(i) ['Var' int2str(i)], s:e, 'UniformOutput', false);
    perf_over_time_chunk = table2array(readall(perf_ds))';
    
    % keep top 10% of performances, which_top tells you the idx of the
    % parameters used
    [perf_top_10(s:e,:), which_top_10(s:e,:)] = maxk(perf_over_time_chunk, round(n/10),2);
end

    % make sure all performances are >=0 by adding the minimum performance
    perf_top_10_pos = perf_top_10 - min(min(perf_top_10,[],'all'),0);

% create empty container for best thetas for each timestep
theta_top_10 = zeros(numel(OF_idx), round(n/10), model.numParams);
% populate with the best sets of theta for each timestep
for t=1:size(theta_top_10,1)
    theta_top_10(t,:,:) = theta_sample(:,which_top_10(t,:))';
end

% calculate cumulative distributions
[cd_top_10, theta_top_10_sorted] = calc_cd_multidim(theta_top_10, perf_top_10_pos);

% calculate information content for each parameter at each timestep
info_content = calc_info_content(cd_top_10, theta_top_10_sorted);

% calculate gradients of the (binned) cumulative distribution
[cd_gradient_breaks, cd_gradient] = calc_cd_gradient(cd_top_10, theta_top_10_sorted);

save(file_log, "cd_gradient", "cd_gradient_breaks", "info_content", "perf_top_10", "-append");