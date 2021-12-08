function time = lambdaC_BCfailure_mean_tta_and_variance(method, k)

n = 6;% number of components

setup_rates
% TODO
lambdaR = (k+1)*0.041665*ones(n,1);
    
setup_case_study;
        
[mtta, timemtta] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol);
[var, timevar] = eval_measure('tta_variance', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol);
                   
results = [lambdaC(1), mtta, var, timemtta, timevar];

fprintf("lambdaC\t\tMTTA\t\tVar\t\tTimeMTTA\tTimeVAR\n");
fprintf("%f\t%f\t%f\t%f\t%f\n", results');

time = timemtta + timevar;
                   
end

