function [cd_gradient,...
          cd_gradient_breaks,...
          info_content] = ...
                          dynia(log_file)

load(log_file);

for j = 1:numel(simdata.OF_idx)
    this_file_name = [o.file_prefix '_' num2str(simdata.OF_idx(j)) '.json'];
end

% create empty containers for top10 performance and parameters at every
% timestep
which_top_10 = zeros(numel(OF_idx), round(n/10));
perf_top_10 = zeros(numel(OF_idx), round(n/10));

% divide this in chunks too otherwise it is too large:
e = 0;
while e < numel(OF_idx)
    s = e + 1; e = min(s + c_size -1, numel(OF_idx));
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

save(file_log, "cd_gradient", "cd_gradient_breaks", "info_content", "-append");