% TODO set parameters   

eta = 0.995

zeta = 0.2

% 0.5 (once every 2 hours)
lambdaR_min = 0.45;
lambdaR_max = 0.55;
lambdaR = lambdaR_min + (lambdaR_max-lambdaR_min)*rand(n,1)

% 0.014 (about once every 72 hours)
lambdaD_min = 0.013;
lambdaD_max = 0.015;
lambdaD = lambdaD_min + (lambdaD_max-lambdaD_min)*rand(n,1)

% 0.002 (once every 500 hours)
lambdaE_min = 0.0015;
lambdaE_max = 0.0025;
lambdaE = lambdaE_min + (lambdaE_max-lambdaE_min)*rand(n,1)

% 0.0014 (about once every month)
lambdaEP_min = 0.0013;
lambdaEP_max = 0.0015;
lambdaEP = lambdaEP_min + (lambdaEP_max-lambdaEP_min)*rand(n,1)