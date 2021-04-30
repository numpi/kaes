function time = lambdaC_BCfailure_moments_tta(method, k)

n = 6;% number of components

%  topology = full((eye(n,n)+sprand(n,n, 1 / n)) > 0)
rng(2)
topology = createTopology(n, 0.2, 'starnoloops');

disp(topology);

lambdaB_min = 18;
lambdaB_max = 20;
lambdaB = lambdaB_min + (lambdaB_max-lambdaB_min)*rand(n,1)

lambdaD_min = 9;
lambdaD_max = 11;
lambdaD = lambdaD_min + (lambdaD_max-lambdaD_min)*rand(n,1)

% lambdaW=342;
lambdaW_min = 320;
lambdaW_max = 350;
lambdaW = lambdaW_min + (lambdaW_max-lambdaW_min)*rand(n,1)

lambdaE_min = 1;
lambdaE_max = 1.5;
lambdaE = lambdaE_min + (lambdaE_max-lambdaE_min)*rand(n,1)

lambdaEP_min = 0.9;
lambdaEP_max = 1.1;
lambdaEP = lambdaEP_min + (lambdaEP_max-lambdaEP_min)*rand(n,1)

pp = symrcm(topology);
disp(topology(pp,pp));

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

shift = 0.01;
it_mult = true;

if strcmp(method, 'ttexpsums2')
	shift = 1e6;
end

lambdaC = (k+1)*0.041665*ones(n,1);
    
% Construct Rs and Ws
[R, W] = BCfailure(n, topology(pp, pp), lambdaB(pp), lambdaC(pp), lambdaD(pp), ...
    lambdaW(pp), lambdaE(pp), lambdaEP(pp));

[m1, timem1] = eval_measure('momentk', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, 'shift', shift, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', 1e-10, 'tol', 1e-4, ...
                       'iterative_mult', it_mult, 'use_sinc', true, ...
                       'interval_report', 10, 'moment', 1);
[m2, timem2] = eval_measure('momentk', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, 'shift', shift, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', 1e-10, 'tol', 1e-4, ...
                       'iterative_mult', it_mult, 'use_sinc', true, ...
                       'interval_report', 10, 'moment', 2);
[m3, timem3] = eval_measure('momentk', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, 'shift', shift, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', 1e-10, 'tol', 1e-4, ...
                       'iterative_mult', it_mult, 'use_sinc', true, ...
                       'interval_report', 10, 'moment', 3);
[m4, timem4] = eval_measure('momentk', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, 'shift', shift, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', 1e-10, 'tol', 1e-4, ...
                       'iterative_mult', it_mult, 'use_sinc', true, ...
                       'interval_report', 10, 'moment', 4);
        

results = [lambdaC(1), m1, m2, m3, m4, timem1, timem2, timem3, timem4];

fprintf("lambdaC\t\tm1\tm2\tm3\tm4\t\tTimeM1\tTimeM2\tTimeM3\tTimeM4\n");
fprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", results');

time = timem4;
                   
end

