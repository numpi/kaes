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

% 0.000114 (once every year)
lambdaE_min = 0.0001135;
lambdaE_max = 0.0001145;
lambdaE = lambdaE_min + (lambdaE_max-lambdaE_min)*rand(n,1)

% 0.015 (little less than once every 3 days)
lambdaEP_min = 0.013;
lambdaEP_max = 0.017;
lambdaEP = lambdaEP_min + (lambdaEP_max-lambdaEP_min)*rand(n,1)