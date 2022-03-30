function of_over_time = calc_of_moving_window(Qsim,Qobs,window,of_name,varargin)

n = size(Qsim,2);
of_over_time = zeros(size(Qsim));

for i = 1:n
    of_over_time(:,i) = matlab.tall.movingWindow(of_name,window,Qobs,Qsim(:,i),varargin{:});
end