function test_examples(N)
%TEST_EXAMPLES 

% Number of test to run
k = 10;

acc_times = [];

for j = 1 : k
    acc_times = [ acc_times, input_n_density01(N, 'ttexpsums2') ];
end

fprintf('Average solution time: %f secs -- variance: \n', mean(acc_times), ...
    var(acc_times));

end

