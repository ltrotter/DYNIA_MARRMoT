n = 10000;
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
of_name = 'of_KGE';
window = 31;

[cd_gradient,...
    info_content,...
    theta_sample,...
    Qsim,...
    perf_over_time] = dynia(m, n, Qobs, of_name, window);

save("test_output.mat","cd_gradient","info_content","theta_sample","Qsim","perf_over_time");
