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

if simdata.n_done < o.n
    % run all the simulations, with the appropriate restarting
    helper_sim_function(file_log, simdata, model, Qobs, o)
end
end

function [] = helper_sim_function(file_log, simdata, model, Qobs, o)

    n_done = simdata.n_done;
    theta_sample = simdata.theta_sample;
    OF_idx = simdata.OF_idx;
    total_time = simdata.total_time;

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

    if o.parallelEval == 1 %parallel evaluation
        poolobj = gcp('nocreate');
        if isempty(poolobj); poolobj = parpool; end
        workers = poolobj.NumWorkers;
        while n_done < o.n
            n_to_do = o.n - n_done;
            if n_to_do < workers * o.chunk_size
                o.chunk_size = ceil(n_to_do/workers);
            end
            spmd_start = tic;
            spmd
%                 this_worker_data = load([o.file_prefix 'workers_data.mat']);
%                 OF_idx = this_worker_data.OF_idx;
%                 model = this_worker_data.model;
%                 numel_OF_idx = this_worker_data.numel_OF_idx;
%                 Qobs = this_worker_data.Qobs;
%                 o = this_worker_data.o;
                this_start_i = n_done + 1 + (labindex-1)*o.chunk_size;
                this_end_i = min(this_start_i + o.chunk_size - 1, size(theta_sample,2));
                this_theta_sample = theta_sample(:,this_start_i:this_end_i);
                this_output = cell(numel(OF_idx),1);
                for j=1:size(this_theta_sample,2)
                    this_theta = this_theta_sample(:,j);
        
                    % run the model with this parameter set
                    model.theta = this_theta; model.run();
                
                    % get the streamflow
                    Qsim = model.get_streamflow();
                     
                    % calculate performance over time
                    perf_over_time = calc_of_moving_window(Qsim, Qobs, o.window, o.step, o.of_name, o.precision_Q+1, o.of_args{:});
            
                    % check the timesteps with performance above the threshold
                    behavioural_steps = perf_over_time > o.perf_thr;
                    OF_idx_behavioural = OF_idx(behavioural_steps);
                    idx_behavioural = find(behavioural_steps);
            
                    % extract performance, fluxes and stores at those timesteps
                    OF_vals = perf_over_time(behavioural_steps);
                    fluxes = calc_avg_at_timesteps(model.fluxes, OF_idx_behavioural, o.window);
                    stores = calc_avg_at_timesteps(model.stores, OF_idx_behavioural, o.window);
        
                    % for each of those steps
                    for i = 1:numel(idx_behavioural)
                        t = idx_behavioural(i);
            
                        % get all the useful data
                        OF_here = round(OF_vals(i), o.precision_OF);
                        fluxes_here = round(fluxes(i,:), o.precision_Q);
                        stores_here = round(stores(i,:), o.precision_Q);
            
                        % create the one line string that will need to be written to
                        % file
                        csv_txt = [num2str(this_theta', '%.9g,'), num2str(OF_here, '%g,'),...
                                   num2str(fluxes_here, '%g,'),  num2str(stores_here, '%g,')];
                        csv_txt(end) = []; csv_txt = [csv_txt, '\n'];
            
                        % append it to what's already existing
                        this_output{t} = [this_output{t}, csv_txt];
                    end
                end
                
                % write outputs of all timesteps to file, after the end of the loop
                for c=1:numel(this_output)
                    if isempty(this_output{c}); continue; end
                    % create the file name
                    file_name = [o.file_prefix '_' num2str(OF_idx(c)) '.csv'];
                    
                    % save the relevant output
                    fileID = fopen(file_name,'a');
                    fprintf(fileID, this_output{c});
                    fclose(fileID);
                end
                
                done_these_chunks = j;
            end
            n_done = n_done + sum([done_these_chunks{:}]);

            % calculate time it took for this chunk
            these_chunks_time = toc(spmd_start);
            total_time = total_time + these_chunks_time;
        
            % save n_done to the log file
            simdata.n_done = n_done;
            simdata.chunk_time(end+1) = these_chunks_time;
            simdata.total_time = total_time;
            save(file_log, "simdata", "-append");
        
            % print to screen if display is active
            if o.display
                msg = ['Simulations done: ', num2str(n_done), '/', num2str(o.n),...
                            ' (', num2str(n_done/o.n*100,'%.2f'), '%) - ',...
                            'Elapsed time = ' num2str(total_time/60,'%.2f'), 'min (+',...
                            num2str(these_chunks_time/60,'%.2f'),')'];
                disp(msg);
            end
        end
        delete(poolobj)

    else % series evaluation
        % before we start the loop, prepare the output cell and start a timer
        this_chunk_output = cell(numel(OF_idx),1);
        chunk_time = tic;
    
        % loop through each set of thetas
        while n_done < o.n
            this_theta = theta_sample(:,n_done + 1);
    
            % run the model with this parameter set
            model.theta = this_theta; model.run();
    
            % get the streamflow
            Qsim = model.get_streamflow();
    
            % calculate performance over time
            perf_over_time = calc_of_moving_window(Qsim, Qobs, o.window, o.step, o.of_name, o.precision_Q+1, o.of_args{:});
    
            % check the timesteps with performance above the threshold
            behavioural_steps = perf_over_time > o.perf_thr;
            OF_idx_behavioural = OF_idx(behavioural_steps);
            idx_behavioural = find(behavioural_steps);
    
            % extract performance, fluxes and stores at those timesteps
            OF_vals = perf_over_time(behavioural_steps);
            fluxes = calc_avg_at_timesteps(model.fluxes, OF_idx_behavioural, o.window);
            stores = calc_avg_at_timesteps(model.stores, OF_idx_behavioural, o.window);
    
            % for each of those steps
            for i = 1:numel(idx_behavioural)
                t = idx_behavioural(i);
    
                % get all the useful data
                OF_here = round(OF_vals(i), o.precision_OF);
                fluxes_here = round(fluxes(i,:), o.precision_Q);
                stores_here = round(stores(i,:), o.precision_Q);
    
                % create the one line string that will need to be written to
                % file
                csv_txt = [num2str(this_theta', '%.9g,'), num2str(OF_here, '%g,'),...
                           num2str(fluxes_here, '%g,'),  num2str(stores_here, '%g,')];
                csv_txt(end) = []; csv_txt = [csv_txt, '\n'];
    
                % append it to what's already existing
                this_chunk_output{t} = [this_chunk_output{t}, csv_txt];
            end
    
            % increase n_done
            n_done = n_done + 1;
    
            % each chunk_size simulations, save the results so far
            if rem(n_done, o.chunk_size) == 0
                for c=1:numel(this_chunk_output)
                    if isempty(this_chunk_output{c}); continue; end
                    % create the file name
                    file_name = [o.file_prefix '_' num2str(OF_idx(c)) '.csv'];
                    
                    % save the relevant output
                    fileID = fopen(file_name,'a');
                    fprintf(fileID, this_chunk_output{c});
                    fclose(fileID);
                end
    
                % empty up this_chunk_output
                this_chunk_output = cell(numel(OF_idx),1);
    
                % calculate time it took for this chunk
                this_chunk_time = toc(chunk_time); chunk_time = tic;
                total_time = total_time + this_chunk_time;
    
                % save n_done to the log file
                simdata.n_done = n_done;
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
                        'theta', [],...
                        'chunk_size', 1000,...
                        'parallelEval', 0);
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
