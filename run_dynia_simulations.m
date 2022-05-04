function [] = run_dynia_simulations(model, Qobs, options)

% get all options
if nargin <3 || isempty(options); options = struct(); end
o = get_dynia_options(options);

% check if this has been started already
file_log   = [o.file_prefix, '.mat'];

% if the log file does not exists or if the user decided to overwrite
if ~isfile(file_log) || o.overwrite

    % create all thetas
    simdata.theta_sample = lhs_sample_par(model,o.n);

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
        % create a csv file for each timestep, this will be used later on.
        for j = 1:numel(OF_idx)
            this_file_name = [o.file_prefix '_' num2str(OF_idx(j)) '.csv'];
            fileID = fopen(this_file_name,'w');
            % write the header braket to each file, this will be helpful later
            fprintf(fileID,'theta, OF, fluxes, stores\n');
            fclose(fileID);
        end
        %simdata.isfirst = ones(size(OF_idx)); % this is used to add commas between entries in the json files
    end
    %isfirst = simdata.isfirst;
    
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
            csv_txt = [mat2str(this_theta), ',', num2str(OF_here, o.precision_OF), ',',...
                       mat2str(fluxes_here', o.precision_Q), ',', mat2str(stores_here', o.precision_Q),'\n'];

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