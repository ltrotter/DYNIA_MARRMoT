function [average_at_timesteps] = calc_avg_at_timesteps(data, timesteps, window)

half_window = floor(window/2);
T = numel(timesteps);

average_at_timesteps = NaN(T,size(data,2));

for i=1:T
    t = timesteps(i);
    average_at_timesteps(i,:) = mean(data(t-half_window:t+half_window,:),1);
end
