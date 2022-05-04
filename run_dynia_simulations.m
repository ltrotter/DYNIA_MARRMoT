function [] = run_dynia_simulations(model, Qobs, options)

% get all options
if nargin <3 || isempty(options); options = struct(); end
o = get_dynia_options(options);

% check if this has been started already
file_log   = [o.file_prefix, '.mat'];

% if the log file does not exists or if the user decided to overwrite
if ~isfile(file_log) || o.overwrite

    % create all thetas
    simdata.theta_sample = lhs_sample_par(model,o.n, o.theta);

    % identify the indices of the OF calculation, based on Qobs
    [~, simdata.OF_idx] = calc_of_moving_window(Qobs,Qobs,o.window,o.step,o.of_name, o.precision_Q+1, o.of_args{:});

    % set that none have happened yet
    simdata.n_done = 0;
    simdata.pc_done = 0;

    % save the options and the thetas to a new log file
    save(file_log, "o", "model", "Qobs", "simdata");

% otherwise, this is a restart
else
    % warn that all options will be loaded (i.e. the ones given are all discarded).
    disp([file_log ' found: options will be loaded.'])
    load(file_log, "o", "model", "Qobs", "simdata");
end

if simdata.n_done < o.n
    % run all the simulations, with the appropriate restarting
    helper_sim_function(file_log, simdata, model, Qobs, o)
end
end

function [] = helper_sim_function(file_log, simdata, model, Qobs, o)

    n_done = simdata.n_done;
    theta_sample = simdata.theta_sample;
    OF_idx = simdata.OF_idx;
    pc_done = simdata.pc_done;

    if n_done == 0
        % create a csv file for each timestep
        % create the header, which is common for all files
        flux_names  = cellfun(@(fn) ['flux_',fn,','], cellstr(model.FluxNames), 'UniformOutput', false);
        store_names = cellfun(@(sn) [sn,','], cellstr(model.StoreNames), 'UniformOutput', false);
        header = [num2str(1:model.numParams, 'theta_%i,'), 'OF_value,',  flux_names{:}, store_names{:}];
        header(end) = []; header = [header, '\n'];
        for j = 1:numel(OF_idx)
            this_file_name = [o.file_prefix '_' num2str(OF_idx(j)) '.csv'];
            fileID = fopen(this_file_name,'w');
            % write the header braket to each file
            fprintf(fileID, header);
            fclose(fileID);
        end
    end

    % loop through each set of thetas
    while n_done < o.n
        this_theta = theta_sample(:,n_done + 1);

        % run the modoel with this parameter set
        model.theta = this_theta; model.run();

        % get the streamflow
        Qsim = model.get_streamflow();

        % calculate performance over time
        perf_over_time = calc_of_moving_window(Qsim, Qobs, o.window, o.step, o.of_name, o.precision_Q+1, o.of_args{:});

        % check the timesteps with performance above the threshold
        behavioural_steps = perf_over_time > o.perf_thr;
        OF_idx_behavioural = OF_idx(behavioural_steps);

        % extract performance, fluxes and stores at those timesteps
        OF_vals = perf_over_time(behavioural_steps);
        fluxes = calc_avg_at_timesteps(model.fluxes, OF_idx_behavioural, o.window);
        stores = calc_avg_at_timesteps(model.stores, OF_idx_behavioural, o.window);

        % for each of those steps
        for t = 1:numel(OF_idx_behavioural)
            % create a cell array containing the useful data
            OF_here = round(OF_vals(t), o.precision_OF);
            fluxes_here = round(fluxes(t,:), o.precision_Q);
            stores_here = round(stores(t,:), o.precision_Q);

            % create the name and text to save as json
            this_OF_idx = OF_idx_behavioural(t);
            file_name = [o.file_prefix '_' num2str(this_OF_idx) '.csv'];
            csv_txt = [num2str(this_theta', '%.9g,'), num2str(OF_here, '%g,'),...
                       num2str(fluxes_here, '%g,'),  num2str(stores_here, '%g,')];
            csv_txt(end) = []; csv_txt = [csv_txt, '\n'];

            % save to the csv file
            fileID = fopen(file_name,'a');
            fprintf(fileID, csv_txt);
            fclose(fileID);
        end

        % increase n_done and add it to the log file
        n_done = n_done + 1;
        simdata.n_done = n_done;

        % display progress at each full percentage points
        if(o.display)
            pc_done_now = n_done/o.n * 100;
            if floor(pc_done_now) > floor(pc_done)
                msg = ['Simulations done: ', num2str(n_done), '/', num2str(o.n),...
                        ' (', num2str(floor(pc_done_now),'%u'), '%).'];
                disp(msg);
            end
        end
        pc_done = pc_done_now;
        simdata.pc_done = pc_done;
        save(file_log, "simdata", "-append");
    end
end

function sample = unif_sample_par(model, n)

    ranges = model.parRanges;
    rd = rand(size(ranges, 1),n);
    sample = rd.*diff(ranges,1,2) + ranges(:,1);

end

function sample = lhs_sample_par(model, n, theta)

    if nargin < 3 || isempty(theta); theta = NaN(model.numParams,1);end
    sample = NaN(model.numParams,n);
    

    ranges = model.parRanges(isnan(theta),:);
    lhs = lhsdesign(n,size(ranges, 1));
    sample(isnan(theta),:) = lhs'.*diff(ranges,1,2) + ranges(:,1);
    sample(~isnan(theta),:) = repmat(theta(~isnan(theta)),1,n);

end

function opts_out = get_dynia_options(opts_in)

    defaultopt = struct('repeats', 2e5, ...
                        'of_name', 'of_KGE', ...
                        'window_size', 31, ...
                        'window_step', 1,...
                        'file_prefix', 'DYNIA', ...
                        'precision_Q', 4, ...
                        'precision_OF', 4,...
                        'perf_thr', 0,...
                        'display', 1,...
                        'overwrite', 0,...
                        'theta', []);
    defaultopt.of_args = cell(0);

    % get options
    opts_out.n = optimget(opts_in, 'repeats', defaultopt, 'fast');
    opts_out.of_name = optimget(opts_in, 'of_name', defaultopt, 'fast');
    opts_out.window = optimget(opts_in, 'window_size', defaultopt, 'fast');
    opts_out.step = optimget(opts_in, 'window_step', defaultopt, 'fast');
    opts_out.file_prefix = optimget(opts_in, 'file_prefix', defaultopt, 'fast');
    opts_out.precision_Q = optimget(opts_in, 'precision_Q', defaultopt, 'fast');
    opts_out.precision_OF = optimget(opts_in, 'precision_OF', defaultopt, 'fast');
    opts_out.of_args = optimget(opts_in, 'of_args', defaultopt, 'fast');
    opts_out.perf_thr = optimget(opts_in, 'perf_thr', defaultopt, 'fast');
    opts_out.display = optimget(opts_in, 'display', defaultopt, 'fast');
    opts_out.overwrite = optimget(opts_in, 'overwrite', defaultopt, 'fast');
    opts_out.theta = optimget(opts_in, 'theta', defaultopt, 'fast');

end

function [of_over_time_non_missing, idx_non_missing] = calc_of_moving_window(Qsim,Qobs,window,step,of_name,precision,varargin)

    if isempty(precision); precision = 3; end
    [T,n] = size(Qsim);

    half_window = floor(window/2);
    idx = (half_window+1):step:(T-half_window-1);
    of_over_time = NaN(numel(idx),n);
    for i = 1:numel(idx)
        rand_perm = rand(2*half_window+1,1)*10^-precision;
        this_idx = idx(i);
        Qobs_in_window = Qobs(this_idx-half_window:this_idx+half_window) + rand_perm;
        if any(Qobs_in_window < 0) || any(isnan(Qobs_in_window)); continue; end
        Qsim_in_window = Qsim(this_idx-half_window:this_idx+half_window,:) + rand_perm;
        for p = 1:n
            of_over_time(i,p) = feval(of_name,Qobs_in_window,Qsim_in_window(:,p),varargin{:});
        end
    end

    missing_timesteps = isnan(of_over_time(:,1));
    of_over_time_non_missing = of_over_time(~missing_timesteps,:);
    idx_non_missing = idx(~missing_timesteps);

end

function [average_at_timesteps] = calc_avg_at_timesteps(data, timesteps, window)

    half_window = floor(window/2);
    T = numel(timesteps);

    average_at_timesteps = NaN(T,size(data,2));

    for i=1:T
        t = timesteps(i);
        average_at_timesteps(i,:) = mean(data(t-half_window:t+half_window,:),1);
    end

end
