function [] = run_dynia_simulation(file_log, simdata, model, Qobs, o)

n_done = simdata.n_done;
theta_sample = simdata.theta_sample;
OF_idx = simdata.OF_idx;
pc_done = simdata.pc_done;

if n_done == 0
    % create a json file for each timestep, this will be used later on.
    for j = 1:numel(OF_idx)
        this_file_name = [o.file_prefix '_' num2str(OF_idx(j)) '.json'];
        fileID = fopen(this_file_name,'w');
        % write the opening braket to each file, this will be helpful later
        fprintf(fileID,'[');
        fclose(fileID);
    end
    simdata.isfirst = ones(size(OF_idx)); % this is used to add commas between entries in the json files
end
isfirst = simdata.isfirst;

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
        % create a struct containing the useful data
        outstruct.theta = this_theta;
        outstruct.OF = round(OF_vals(t), o.precision_OF);
        outstruct.fluxes = round(fluxes(t,:), o.precision_Q);
        outstruct.stores = round(stores(t,:), o.precision_Q);

        % create the name and text to save as json
        this_OF_idx = OF_idx_behavioural(t);
        file_name = [o.file_prefix '_' num2str(this_OF_idx) '.json'];
        jsontxt = jsonencode(outstruct);

        % save to the json file
        fileID = fopen(file_name,'a');
        % this if-else is because we need to add commas between each entry
        if isfirst(OF_idx == this_OF_idx)
            isfirst(OF_idx == this_OF_idx) = 0;
        else
            fprintf(fileID,',');
        end
        fprintf(fileID,jsontxt);
        fclose(fileID);
    end

    % increase n_done and add it to the log file
    n_done = n_done + 1;
    simdata.n_done = n_done;
    simdata.isfirst = isfirst;

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

% close the bracket for all files at the end of the operation
for j = 1:numel(OF_idx)
    this_file_name = [o.file_prefix '_' num2str(OF_idx(j)) '.json'];
    fileID = fopen(this_file_name,'a');
    % write the opening braket to each file, this will be helpful later
    fprintf(fileID,']');
    fclose(fileID);
end