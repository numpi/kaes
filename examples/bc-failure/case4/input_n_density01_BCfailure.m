function time = input_n_density01_BCfailure(n, method, topology)
%
% Constructs the example decribed in the paper "New Method for 
% the Evaluation of (Conditional) Time To Failure and its Moments
% Trough Implicit Reward Structure", by G. Masetti, L. Robol, S. Chiaradonna,
%   F. Di Giandomenico, submitted.

%  topology = full((eye(n,n)+sprand(n,n, 1 / n)) > 0)
if ~exist('topology', 'var')
    topology = createTopology(n, 0.2, 'starnoloops');
end

rng(2)

disp(topology);

setup_rates;

pp = symrcm(topology);
disp(topology(pp,pp));
  
% Construct Rs and Ws
[R, W] = BCfailure(n, topology, eta, zeta, lambdaR, lambdaD, lambdaE, lambdaEP, l);

P = eye(5); P(1:3,1:3) = 0;
w = cell(1, n);
for j = 1 : n
	w{j} = tt_matrix(P);
end
absorbing_states = ktt_kron(w{:});

pi0 = ktt_ej(5*ones(1,n), ones(1,n));

v = [ 0 ; 0 ; 0 ; 1 ; 1 ];
w = cell(1, n); for j = 1 : n; w{j} = tt_tensor(v); end

r = round(ktt_ones(5*ones(1,n)) - ktt_kron(w{:}), 1e-8);

% Compute the measure
% if n <= 6
% 	m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
% 				 'algorithm', 'spantree', ...
% 				 'absorbing_states', absorbing_states);
% end

shift = 0.01;
it_mult = true;

if strcmp(method, 'ttexpsums2')
	shift = 1e6;
end

[m, time] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
					   'algorithm', method, 'shift', shift, ...
					   'absorbing_states', absorbing_states, ...
					   'ttol', 1e-10, 'tol', 1e-4, ...
					   'iterative_mult', it_mult, 'use_sinc', true, ...
					   'interval_report', 10);

end
