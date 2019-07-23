function time = input_n_mttf(n, method, topology, it_mult)
%INPUT_N_MTTF 

if ~exist('it_mult', 'var')
	it_mult = false;
end

%  topology = full((eye(n,n)+sprand(n,n, 1 / n)) > 0)
if ~exist('topology', 'var')
    topology = createTopology(n, 0.2, 'starnoloops');
end

lambdah = 0.1 * (1 : n);
lambdahp = (1 : n);
lambdas = (1 : n);

pp = symrcm(topology);

[R, W, ~] = mttfCaseStudy(n, topology(pp, pp), lambdah(pp), ...
	lambdahp(pp), lambdas(pp));

absorbing_states = 3 * ones(1, n);

pi0 = ktt_ej(3*ones(1,n), ones(1,n));
r   = round(ktt_ones(3*ones(1,n)) - ...
			ktt_ej(3*ones(1,n),3*ones(1,n)), 1e-8);

% Compute the measure
if n <= 10
	m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
				 'algorithm', 'spantree', ...
				 'absorbing_states', absorbing_states);
end

shift = 0;

if strcmp(method, 'ttexpsums2')
	shift = 3e1;
end

[m, time] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
					   'algorithm', method, 'shift', shift, ...
					   'absorbing_states', absorbing_states, ...
					   'ttol', 1e-10, 'tol', 1e-3, ...
					   'iterative_mult', it_mult, 'use_sinc', false);
				   
% pause
% [m, time] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
%    'algorithm', 'ttexpsumst', 'shift', 0, ...
%    'absorbing_states', absorbing_states, ...
%    'ttol', 1e-8, 'tol', 1e-3, 'iterative_mult', it_mult);

end

