function [] = run_dynia_simulations(model, Qobs, options)

% get all options
if nargin <3 || isempty(options); options = struct(); end
o = get_dynia_options(options);
if isstruct(o.windows)
    step  = o.windows.step;
    width = o.windows.size;

    n_of_steps = length(Qobs);
    w_start = 1:step:(n_of_steps-width);
    w_end   = width:step:n_of_steps; 

    windows = [w_start', w_end'];
    o.windows = windows;
end

% check if this has been started already
file_log   = [o.file_prefix, '.mat'];

% if the log file does not exists or if the user decided to overwrite
if ~isfile(file_log) || o.overwrite

    % create all thetas
    simdata.theta_sample = lhs_sample_par(model,o.n, o.theta, 1234);

    % identify the indices of the OF calculation, based on Qobs
    [~, nonmissingwindows] = calc_of_moving_window(Qobs,Qobs,o.windows,o.of_name, o.precision_Q+1, o.of_args{:});
    simdata.OF_idx = o.windows(nonmissingwindows, :);


    % set that none have happened yet
    simdata.to_do = ones(1,o.n);
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

if any(simdata.to_do)
    % prepare output (i.e. list of 10% of best performing sets for each ts)
    n_to_keep       = ceil(o.n * o.pc_top);
    top_performance = o.sign * inf(n_to_keep, size(simdata.OF_idx, 1));
    if all(simdata.to_do)
        % create a new mat file for each timestep to store the top 10% of
        % simulation results
        output = cell(size(top_performance, 1), 1);
        for j = 1:size(simdata.OF_idx, 1)
            this_file_name = [o.file_prefix, '_',...
                              num2str(simdata.OF_idx(j, 1)), '-',...
                              num2str(simdata.OF_idx(j, 2)), '.mat'];
            save(this_file_name, 'output')
        end
    end
    while any(simdata.to_do)
        % run a single chunk
        chunk_time = tic;
        if o.parallelEval == 0
        [top_perf, output, theta_done] = ...
            run_chunk_simulations(simdata, top_performance, model, Qobs, o);
        else
            poolobj = gcp('nocreate');
            if isempty(poolobj); parpool; end
            %workers = poolobj.NumWorkers;
            ids_to_do = find(simdata.to_do);
            spmd
                this_start_i = 1 + (labindex-1)*o.chunk_size;
                this_end_i = min(this_start_i + o.chunk_size - 1, size(ids_to_do,2));
                this_theta_ids = ids_to_do(this_start_i:this_end_i);
                this_simdata = simdata;
                this_simdata.to_do = zeros(size(simdata.to_do));
                this_simdata.to_do(this_theta_ids) = 1;

                [top_perf_par, output_par, theta_done_par] = ...
                    run_chunk_simulations(this_simdata, top_performance, model, Qobs, o);
            end

            [top_perf, output, theta_done] = ...
                combine_from_workers(top_perf_par, output_par, theta_done_par,o.sign);

        end
        simdata.to_do(theta_done) = 0;
        top_performance = top_perf.perf;
        save_chunk_data(output, top_perf, o, simdata.OF_idx)
        simdata.chunk_time(end+1) = toc(chunk_time);
        simdata.total_time = sum(simdata.chunk_time);
        if o.display; print_summary(simdata, o); end
    end

    % create the header, which is common for all files
    flux_names  = cellfun(@(fn) ['flux_',fn,','], cellstr(model.FluxNames), 'UniformOutput', false);
    store_names = cellfun(@(sn) [sn,','], cellstr(model.StoreNames), 'UniformOutput', false);
    header = ['OF_value,',  num2str(1:model.numParams, 'theta_%i,'),flux_names{:}, store_names{:}];
    header(end) = []; header = [header, '\n'];
    for i=1:size(simdata.OF_idx,1)
        this_file = [o.file_prefix, '_', num2str(simdata.OF_idx(i, 1)),...
                     '-', num2str(simdata.OF_idx(i, 2))];
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

function [top_perf, output, theta_done] = ...
    run_chunk_simulations(simdata, top_performance, model, Qobs, o)

    theta_sample = simdata.theta_sample;
    OF_idx = simdata.OF_idx;
    theta_done = [];

    % before we start the loop, prepare the output cell
    this_chunk_output = cell(o.chunk_size,size(OF_idx, 1));

    % set that this is a new chunk
    chunk_id = 0;
    % and all top performance values belong to the previous chunks.
    where_top_performance = zeros(size(top_performance));

    % loop through each set of thetas
    while chunk_id < o.chunk_size && any(simdata.to_do); chunk_id = chunk_id +1;
        % get the first undone parameter set
        this_theta_id = find(simdata.to_do, 1, 'first');
        this_theta = theta_sample(:,this_theta_id);

        % run the model with this parameter set
        model.theta = this_theta; model.run();

        % get the streamflow
        Qsim = model.get_streamflow();

        % calculate performance over time
        perf_over_time = calc_of_moving_window(Qsim, Qobs, o.windows, o.of_name, o.precision_Q+1, o.of_args{:});

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
        OF_idx_to_keep = OF_idx(steps_to_keep, :);
        fluxes = calc_avg_in_windows(model.fluxes, OF_idx_to_keep);
        stores = calc_avg_in_windows(model.stores, OF_idx_to_keep);

        tmp = find(steps_to_keep);
        % for each of those steps
        for i = 1:size(OF_idx_to_keep, 1)
            
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
        simdata.to_do(this_theta_id) = 0;
        theta_done(end+1) = this_theta_id;
    end

    % format output
    top_perf.perf = top_performance;
    top_perf.idx = where_top_performance;
    output = this_chunk_output;
end

function [] = save_chunk_data(chunk_output, top_perf, o, OF_idx)

    where_top_performance = top_perf.idx;

    % for each timestep
    for i=1:size(chunk_output, 2)
        if(all(where_top_performance(:,i)==0));continue;end

        % load the existing results
        this_t_filename = [o.file_prefix, '_', num2str(OF_idx(i, 1)), ...
                           '-', num2str(OF_idx(i, 2))];
        load(this_t_filename, 'output');
        old_output = output;
        
        % build the new output
        output_new = old_output; % start with the existing values
        to_change = where_top_performance(:,i) ~= 0; % where the top performance values are in this chunk
        to_change_with = where_top_performance(to_change, i); % find the output string to substitute in
        output_new(to_change) = chunk_output(to_change_with, i); % do the substitution

        % save to the file
        output = output_new;
        save(this_t_filename, 'output');
    end
end

function [] = print_summary(simdata, o)
    
    n_done = sum(simdata.to_do == 0);
    total_time = simdata.total_time;
    last_chunk_time = simdata.chunk_time(end);
    msg = ['Simulations done: ', num2str(n_done), '/', num2str(o.n),...
                ' (', num2str(n_done/o.n*100,'%.2f'), '%) - ',...
                'Elapsed time = ' num2str(total_time/60,'%.2f'), 'min (+',...
                num2str(last_chunk_time/60,'%.2f'),')'];
    disp(msg);

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
    defaultopt.windows.size = 31;
    defaultopt.windows.step = 1;

    % get options
    opts_out.n = optimget(opts_in, 'repeats', defaultopt, 'fast');
    opts_out.of_name = optimget(opts_in, 'of_name', defaultopt, 'fast');
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

    opts_out.windows = optimget(opts_in, 'windows', defaultopt, 'fast');

end

function [of_over_time_non_missing, idx_non_missing] = calc_of_moving_window(Qsim,Qobs,windows,of_name,precision,varargin)

    if isempty(precision); precision = 3; end

    n = size(Qsim, 2);  % number of parameter sets (this is actually always 1, it used to be different)
    t = size(windows,1); % number of windows

    of_over_time = NaN(t,n);
    for i = 1:t
        this_window = windows(i, :);
        Qobs_in_window_raw = Qobs(this_window(1) : this_window(2));
        if any(Qobs_in_window_raw < 0) || any(isnan(Qobs_in_window_raw)); continue; end
        this_window_size = length(Qobs_in_window_raw);
        rand_perm = rand(this_window_size,1)*10^-precision; 
        Qobs_in_window = Qobs_in_window_raw + rand_perm;
        Qsim_in_window = Qsim(this_window(1) : this_window(2), :) + rand_perm;
        for p = 1:n
            of_over_time(i,p) = feval(of_name,Qobs_in_window,Qsim_in_window(:,p),varargin{:});
        end
    end

    missing_timesteps = isnan(of_over_time(:,1));
    of_over_time_non_missing = of_over_time(~missing_timesteps,:);
    idx_non_missing = find(~missing_timesteps);

end

function [average_in_windows] = calc_avg_in_windows(data, windows)

    t = size(windows,1);

    average_in_windows = NaN(t,size(data,2));

    for i=1:t
        average_in_windows(i,:) = mean(data(windows(i,1):windows(i,2),:),1);
    end

end

function  [top_perf, output, theta_done] = ...
                combine_from_workers(top_perf_par, output_par, theta_done_par, perf_sign)
    
    % get everything from the first worker, this will later be update with
    % data from the other workers
    workers = numel(top_perf_par);
    top_perf = top_perf_par{1};
    top_performance = top_perf.perf;
    where_top_performance = top_perf.idx;
    output = output_par{1};
    theta_done = theta_done_par{1};

    for w=2:workers

        this_w_perf = top_perf_par{w};

        this_w_performance = this_w_perf.perf;
        where_this_w_performance = this_w_perf.idx;

        % end here if all values are from previous chunks
        if all(where_this_w_performance == 0); continue; end 

        % all values that come from previous chunks can be set to really bad
        this_w_performance(where_this_w_performance == 0) = perf_sign * inf;

        % arrange the values in descending order
        [~, order] = sort(perf_sign * this_w_performance, 1);
        real_order = sub2ind(size(order), order, repmat(1:size(order,2), size(order,1), 1));
        this_w_performance = this_w_performance(real_order);
        where_this_w_performance = where_this_w_performance(real_order);

        % remove rows made up completely of inf
        allinf = all(isinf(this_w_performance),2);
        this_w_performance(allinf,:) = [];
        where_this_w_performance(allinf,:) = [];

        % get the output from this worker
        this_w_output = output_par{w};

        % add it to the output from the previous workers
        n_output = size(output,1);
        output = [output; this_w_output];

        % update the list of parameters done
        theta_done = [theta_done, theta_done_par{w}];

        for s=1:size(this_w_performance,1)
            this_performance = this_w_performance(s,:);
            where_this_performance = where_this_w_performance(s,:);

            % check the steps where the new performance is better than the
            % worst one in the top 10%
            [worst_top_performance, subs] = (max(top_performance*perf_sign, [], 1));
            steps_to_keep  = this_performance * perf_sign < worst_top_performance & where_this_performance ~= 0;
            % end here if you've got nothing to keep
            if all(~steps_to_keep); continue; end

            % this is all needed to do the substitutions into
            % 'top_performance'
            [a,b] = size(top_performance);
            subs_real_tmp = (((a * (1:b)) - a) + subs);
            subs_real = subs_real_tmp(steps_to_keep);

            % update the top performance for the next step, so that we
            % don't do unnecessary calculations in case this one is better
            % than the next theta
            this_OF_vals = this_performance(steps_to_keep);
            top_performance(subs_real) = this_OF_vals;

            %this tells us where the value of top_performance is found (0=previous chunk)
            where_top_performance(subs_real) = where_this_performance(steps_to_keep) + n_output;
        end     
    end

    % put everything together in top_perf before ending the function
    top_perf.perf = top_performance;
    top_perf.idx = where_top_performance;
end