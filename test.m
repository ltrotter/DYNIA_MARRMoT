m = m_07_gr4j_4p_2s();
load MARRMoT_example_data.mat
input_climatology.precip   = data_MARRMoT_examples.precipitation;                   % Daily data: P rate  [mm/d]
input_climatology.temp     = data_MARRMoT_examples.temperature;                     % Daily data: mean T  [degree C]
input_climatology.pet      = data_MARRMoT_examples.potential_evapotranspiration;    % Daily data: Ep rate [mm/d]
input_climatology.delta_t  = 1;

m.S0 = zeros(m.numStores,1);
m.solver_opts = m.default_solver_opts();
m.input_climate = input_climatology;

Qobs = data_MARRMoT_examples.streamflow;

opts.repeats = 250;
opts.chunk_size = 25;
opts.window_size = 31;
opts.window_step = 7;
opts.of_name = 'of_bias_penalised_log';
opts.of_args = {[],'of_mean_hilo_root5_KGE'};
opts.file_prefix = 'test_dynia';

dynia(m, Qobs, opts);
