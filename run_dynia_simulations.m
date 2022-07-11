function [] = run_dynia_simulations(model, Qobs, options)

% get all options
if nargin <3 || isempty(options); options = struct(); end
o = get_dynia_options(options);

% check if this has been started already
file_log   = [o.file_prefix, '.mat'];

% if the log file does not exists or if the user decided to overwrite
if ~isfile(file_log) || o.overwrite

    % create all thetas
    simdata.theta_sample = lhs_sample_par(model,o.n, o.theta, 1234);

    % identify the indices of the OF calculation, based on Qobs
    [~, simdata.OF_idx] = calc_of_moving_window(Qobs,Qobs,o.window,o.step,o.of_name, o.precision_Q+1, o.of_args{:});

    % set that none have happened yet
    %simdata.n_done = 0;
    simdata.done = zeros(1,o.n);
    simdata.total_time = 0;
    simdata.chunk_time = [];

    % save the options and the thetas to a new log file
    save(file_log, "o", "model", "Qobs", "simdata");

% otherwise, this is a restart
else
    % warn that all options will be loaded (i.e. the ones given are all discarded).
    disp([file_log ' found: options will be loaded.'])
    load(file_log, "o", "model", "Qobs", "simdata");
end

if any(~simdata.done)
    if o.parallelEval == 0
        % run all the simulations, with the appropriate restarting
        helper_sim_function(file_log, simdata, model, Qobs, o)
    else 

    % create the header, which is common for all files
    flux_names  = cellfun(@(fn) ['flux_',fn,','], cellstr(model.FluxNames), 'UniformOutput', false);
    store_names = cellfun(@(sn) [sn,','], cellstr(model.StoreNames), 'UniformOutput', false);
    header = ['OF_value,',  num2str(1:model.numParams, 'theta_%i,'),flux_names{:}, store_names{:}];
    header(end) = []; header = [header, '\n'];
    for i=1:numel(simdata.OF_idx)
        t = simdata.OF_idx(i);
        this_file = [o.file_prefix, '_', num2str(t)];
        load(this_file, 'output');

        csv_file = [this_file, '.csv'];
        fileID = fopen(csv_file,'w');
        % write the header to each file
        fprintf(fileID, header);

        % write the output under the header for each file
        cellfun(@(out) fprintf(fileID, [out, '\n']), output);
        fclose(fileID);

        % delete the .mat file
        delete([this_file, '.mat']);
    end

end
end

function [] = helper_sim_function(file_log, simdata, model, Qobs, o)

    theta_sample = simdata.theta_sample;
    OF_idx = simdata.OF_idx;
    total_time = simdata.total_time;
    
    % prepare output (i.e. list of 10% of best performing sets for each ts)
    n_to_keep       = ceil(o.n * o.pc_top);
    top_performance = o.sign * inf(n_to_keep, numel(simdata.OF_idx));
    
    all_performance = zeros(o.n, numel(OF_idx));

    if all(~simdata.done)
        % create a new mat file for each timestep to store the top 10% of
        % simulation results
        output = cell(n_to_keep, 1);
        for j = 1:numel(OF_idx)
            this_file_name = [o.file_prefix '_' num2str(OF_idx(j)) '.mat'];
            save(this_file_name, 'output')
        end
    end

    % before we start the loop, prepare the output cell and start a timer
    this_chunk_output = cell(o.chunk_size,numel(OF_idx));
    chunk_time = tic;

    % set that this is a new chunk
    chunk_id = 0;
    % and all top performance values belong to the previous chunks.
    where_top_performance = zeros(size(top_performance));

    % loop through each set of thetas
    while any(~simdata.done); chunk_id = chunk_id +1;
        % get the first undone parameter set
        this_theta_id = find(~simdata.done, 1, 'first');
        this_theta = theta_sample(:,this_theta_id);

        % run the model with this parameter set
        model.theta = this_theta; model.run();

        % get the streamflow
        Qsim = model.get_streamflow();

        % calculate performance over time
        perf_over_time = calc_of_moving_window(Qsim, Qobs, o.window, o.step, o.of_name, o.precision_Q+1, o.of_args{:});
        all_performance(this_theta_id, :) = perf_over_time;

        % check the steps where the new performance is better than the
        % worst one in the top 10%
        [worst_top_performance, subs] = (max(top_performance*o.sign, [], 1));
        steps_to_keep  = perf_over_time * o.sign < worst_top_performance';

        % this is all needed to do the substitutions into
        % 'top_performance'
        [a,b] = size(top_performance);
        subs_real_tmp = (((a * (1:b)) - a) + subs);
        subs_real = subs_real_tmp(steps_to_keep);

        % update the top performance for the next step, so that we
        % don't do unnecessary calculations in case this one is better
        % than the next theta
        OF_vals = perf_over_time(steps_to_keep);
        top_performance(subs_real) = OF_vals;
        where_top_performance(subs_real) = chunk_id; %this tells us where the value of top_performance is found (0=previous chunk)

        % get values of stores and fluxes foe each of the ones that
        % needs keeping
        OF_idx_to_keep = OF_idx(steps_to_keep);
        fluxes = calc_avg_at_timesteps(model.fluxes, OF_idx_to_keep, o.window);
        stores = calc_avg_at_timesteps(model.stores, OF_idx_to_keep, o.window);

        % for each of those steps
        for i = 1:numel(OF_idx_to_keep)
            tmp = find(steps_to_keep);
            t = tmp(i);

            % get all the useful data
            OF_here = round(OF_vals(i), o.precision_OF);
            fluxes_here = round(fluxes(i,:), o.precision_Q);
            stores_here = round(stores(i,:), o.precision_Q);

            % create the one line string that will need to be written to
            % file
            csv_txt = [num2str(OF_here, '%g,'),...
                       num2str(this_theta', '%.9g,'),...
                       num2str(fluxes_here, '%g,'),  num2str(stores_here, '%g,')];
            csv_txt(end) = [];% csv_txt = [csv_txt, '\n'];

            % append it to what's already existing
            this_chunk_output{chunk_id, t} = csv_txt;
        end

        % update the list of done sets
        simdata.done(this_theta_id) = 1;
        n_done = sum(simdata.done);

        % each chunk_size simulations, save the results so far
        if chunk_id == o.chunk_size || all(simdata.done)
            % for each timestep
            for i=1:size(this_chunk_output, 2)
                if(all(where_top_performance(:,i)==0));continue;end
                t = OF_idx(i); %get the timestep number

                % load the existing results
                this_t_filename = [o.file_prefix, '_', num2str(t)];
                load(this_t_filename, 'output');
                old_output = output;
                
                % build the new output
                output_new = old_output; % start with the existing values
                to_change = where_top_performance(:,i) ~= 0; % where the top performance values are in this chunk
                to_change_with = where_top_performance(to_change, i); % find the output string to substitute in
                output_new(to_change) = this_chunk_output(to_change_with, i); % do the substitution

                % save to the file
                output = output_new;
                save(this_t_filename, 'output');
            end

            % empty up this_chunk_output
            this_chunk_output = cell(o.chunk_size,numel(OF_idx));

            % set that this is a new chunk
            chunk_id = 0;
            % and all top performance values belong to the previous chunks.
            where_top_performance = zeros(size(top_performance));

            % calculate time it took for this chunk
            this_chunk_time = toc(chunk_time); chunk_time = tic;
            total_time = total_time + this_chunk_time;

            % save n_done to the log file
            simdata.chunk_time(end+1) = this_chunk_time;
            simdata.total_time = total_time;
            save(file_log, "simdata", "-append");

            % print to screen if display is active
            if o.display
                msg = ['Simulations done: ', num2str(n_done), '/', num2str(o.n),...
                            ' (', num2str(n_done/o.n*100,'%.2f'), '%) - ',...
                            'Elapsed time = ' num2str(total_time/60,'%.2f'), 'min (+',...
                            num2str(this_chunk_time/60,'%.2f'),')'];
                disp(msg);
            end
        end
    end
end

function sample = unif_sample_par(model, n)

    ranges = model.parRanges;
    rd = rand(size(ranges, 1),n);
    sample = rd.*diff(ranges,1,2) + ranges(:,1);

end

function sample = lhs_sample_par(model, n, theta, seed)

    if nargin == 4 || ~isempty(seed); rng(seed); end

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
                        'theta', [],...
                        'chunk_size', 1000,...
                        'parallelEval', 0, ...
                        'OF_sign',-1,...
                        'pc_top_performance', .1);
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
    opts_out.chunk_size = optimget(opts_in, 'chunk_size', defaultopt, 'fast');
    opts_out.parallelEval = optimget(opts_in, 'parallelEval', defaultopt, 'fast');
    opts_out.sign = optimget(opts_in, 'OF_sign', defaultopt, 'fast');
    opts_out.pc_top = optimget(opts_in, 'pc_top_performance', defaultopt, 'fast');

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

function [] = save_tmp_data(tmp_file, done_ids, output)
    save(tmp_file, "done_ids", "output");
end