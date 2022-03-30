function of_over_time = calc_of_moving_window(Qsim,Qobs,window,of_name,varargin)

[T,n] = size(Qsim);
of_over_time = NaN(T,n);

step = floor(window/2);
for i = (step+1):(T-step-1)
    Qobs_in_window = Qobs(i-step:i+step);
    if any(Qobs_in_window < 0) || any(isnan(Qobs_in_window)); continue; end
    Qsim_in_window = Qsim(i-step:i+step,:);
    for p = 1:n
        of_over_time(i,p) = feval(of_name,Qobs_in_window,Qsim_in_window(:,p),varargin{:});
    end
end