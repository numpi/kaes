function time = varying_pC_BCfailure_mean_tta_and_variance(method)

n = 3;% number of components

nchanges = 5;% number of different value of pC between 0 and 1-pEP

%  topology = full((eye(n,n)+sprand(n,n, 1 / n)) > 0)
topology = createTopology(n, 0.2, 'starnoloops');

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

% Compute the measure
pEP = 0.2222;

results = zeros(nchanges+1,3);

for k = 0 : nchanges
    pC = k/nchanges*(1-pEP);
    pD = 1-pC-pEP;
    % Construct Rs and Ws
    [R, W] = BCfailure(n, topology(pp, pp), lambdaB(pp), lambdaC(pp), lambdaD(pp), ...
        lambdaW(pp), lambdaE(pp), lambdaEP(pp), pC, pD, pEP);
    [mtta, time] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
                           'algorithm', method, 'shift', shift, ...
                           'absorbing_states', absorbing_states, ...
                           'ttol', 1e-10, 'tol', 1e-4, ...
                           'iterative_mult', it_mult, 'use_sinc', true, ...
                           'interval_report', 10);
    [var, time] = eval_measure('tta_variance', pi0, r, R, W, 'debug', true, ...
                           'algorithm', method, 'shift', shift, ...
                           'absorbing_states', absorbing_states, ...
                           'ttol', 1e-10, 'tol', 1e-4, ...
                           'iterative_mult', it_mult, 'use_sinc', true, ...
                           'interval_report', 10);
    results(k+1,:) = [pC, mtta, var];
end

fprintf("pC\tMTTF\tVar\n");
fprintf("%f\t%f\t%f\n", results');
                   
end

