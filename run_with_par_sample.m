function Qout = run_with_par_sample(m, theta_sample)

n = size(theta_sample,2);
Qout = zeros(size(m.input_climate,1), n);

for i = 1:n
    m.theta = theta_sample(:,i);
    Qout(:,i) = m.get_streamflow();
end