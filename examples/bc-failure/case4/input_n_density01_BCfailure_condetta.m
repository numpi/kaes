function time = input_n_density01_BCfailure_condetta(n, method, casenumber)
%
% Constructs the example decribed in the paper "New Method for 
% the Evaluation of (Conditional) Time To Failure and its Moments
% Trough Implicit Reward Structure", by G. Masetti, L. Robol, S. Chiaradonna,
%   F. Di Giandomenico, submitted.

setup_rates;
setup_case_study;

switch casenumber
    case 1
        % Catastrophic case
        conditional_indices = 4 * ones(1, n);
    case 2
        % All in the benign case
        conditional_indices = [ 4 * ones(1,l), 5 * ones(1, n-l) ];
    case 3
        % Catastrophic + exactly one benign case, for i > l
        conditional_indices = zeros(0, n);
        for i = l+1 : n
            w = 4 * ones(1, n); 
            w(i) = 5;
            conditional_indices = [ conditional_indices ; w ];
        end
    case 4
        % Catastrophic + exactly two benign cases, also for i > l
        conditional_indices = zeros(0, n);
        for i1 = l+1 : n
            for i2 = i1 + 1 : n
                w = 4 * ones(1, n); 
                w(i1) = 5;
                w(i2) = 5;
                conditional_indices = [ conditional_indices ; w ];
            end
        end
end

[m, time] = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
					   'algorithm', method, 'batch_size', 2, ...
					   'absorbing_states', absorbing_states, ...
					   'conditional_indices', conditional_indices, ...
					   'ttol', ttol, 'tol', tol);

end

