function time = input_n_density01_BCfailure_condetta(n, method, topology)
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

lambdaB_min = 100;
lambdaB_max = 200;
lambdaB = lambdaB_min + (lambdaB_max-lambdaB_min)*rand(n,1)

lambdaC_min = 0.1;
lambdaC_max = 0.2;
lambdaC = lambdaC_min + (lambdaC_max-lambdaC_min)*rand(n,1)

lambdaD_min = 1;
lambdaD_max = 1.5;
lambdaD = lambdaD_min + (lambdaD_max-lambdaD_min)*rand(n,1)

lambdaW_min = 1;
lambdaW_max = 1.5;
lambdaW = lambdaW_min + (lambdaW_max-lambdaW_min)*rand(n,1)

lambdaE_min = 1;
lambdaE_max = 1.5;
lambdaE = lambdaE_min + (lambdaE_max-lambdaE_min)*rand(n,1)

lambdaEP_min = 1;
lambdaEP_max = 1.5;
lambdaEP = lambdaEP_min + (lambdaEP_max-lambdaEP_min)*rand(n,1)

  
% Construct Rs and Ws
[R, W] = BCfailure(n, topology, lambdaB, lambdaC, lambdaD, ...
	lambdaW, lambdaE, lambdaEP);

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

conditional_indices = 4 * ones(1, n);

% Compute the measure
if n <= 6
	m = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
				 'algorithm', 'spantree', ...
				 'conditional_indices', conditional_indices, ...
				 'absorbing_states', absorbing_states);
end

shift = 0;
it_mult = false;

if strcmp(method, 'ttexpsums2')
	shift = 2e6;
end

[m, time] = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
					   'algorithm', method, 'shift', shift, ...
					   'absorbing_states', absorbing_states, ...
					   'conditional_indices', conditional_indices, ...
					   'ttol', 1e-10, 'tol', 1e-4, ...
					   'iterative_mult', it_mult, 'use_sinc', false, ...
					   'interval_report', 10);

end

