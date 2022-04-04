function [idx_non_missing, of_over_time_non_missing] = calc_of_moving_window(Qsim,Qobs,window,step,of_name,varargin)

[T,n] = size(Qsim);

half_window = floor(window/2);
idx = (half_window+1):step:(T-half_window-1);
of_over_time = NaN(numel(idx),n);
for i = 1:numel(idx)
    this_idx = idx(i);
    Qobs_in_window = Qobs(this_idx-half_window:this_idx+half_window);
    if any(Qobs_in_window < 0) || any(isnan(Qobs_in_window)); continue; end
    Qsim_in_window = Qsim(this_idx-half_window:this_idx+half_window,:);
    for p = 1:n
        of_over_time(i,p) = feval(of_name,Qobs_in_window,Qsim_in_window(:,p),varargin{:});
    end
end

missing_timesteps = isnan(of_over_time(:,1));
of_over_time_non_missing = of_over_time(~missing_timesteps,:);
idx_non_missing = idx(~missing_timesteps);