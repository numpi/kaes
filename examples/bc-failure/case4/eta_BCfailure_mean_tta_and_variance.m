function time = eta_BCfailure_mean_tta_and_variance(method, k)

n = 3;% number of components

setup_rates

eta = (k+1)*0.001 + 0.989;
    
setup_case_study;
        
[mtta, timemtta] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol);
[var, timevar] = eval_measure('tta_variance', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol);
                   
results = [eta, mtta, var, timemtta, timevar];

fprintf("eta\t\tMTTA\t\tVar\t\tTimeMTTA\tTimeVAR\n");
fprintf("%f\t%f\t%f\t%f\t%f\n", results');

time = timemtta + timevar;
                   
end

