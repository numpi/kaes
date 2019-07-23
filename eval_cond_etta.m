function [m, time] = eval_cond_etta(pi0, R, W, absorbing_states, ...
				algorithm, debug, tol, ttol, shift, iterative_mult, ...
				use_sinc, interval_report)
			
k = length(R);
n = arrayfun(@(i) size(R{i}, 1), 1 : k);
			
Q = ktt_infgen(R, W);
			
% Compute the column vector a with the exit rates for the first 
% absorbing state
a = Q * ktt_ej(n, absorbing_states(1, :));

% Compute the inverse function Qh \ a
[m1, time1, v] = eval_inv(pi0, a, R, W, absorbing_states, algorithm, debug, ...
						 tol, ttol, shift, iterative_mult, use_sinc, ...
						 interval_report);
[m2, time2, ~] = eval_inv(pi0, -v, R, W, absorbing_states, algorithm, debug, ...
						 tol, ttol, shift, iterative_mult, use_sinc, ...
						 interval_report);

m = m2 / m1;

time = time1 + time2;
		 
			
end