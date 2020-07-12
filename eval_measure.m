function [m, time] = eval_measure(fun, pi0, rewards, R, W, varargin)
%EVAL_MEASURE 

p = inputParser;

addOptional(p, 'algorithm', 'auto');
addOptional(p, 'absorbing_states', []);
addOptional(p, 'debug', true);
addOptional(p, 'tol', 1e-4);
addOptional(p, 'ttol', 1e-8);
addOptional(p, 'shift', 0);
addOptional(p, 'iterative_mult', false);
addOptional(p, 'use_sinc', false);
addOptional(p, 'interval_report', 50);
addOptional(p, 'conditional_indices', []);
addOptional(p, 'x0', []);
addOptional(p, 'anderson', false);

parse(p, varargin{:});

algorithm = p.Results.algorithm;
absorbing_states = p.Results.absorbing_states;
debug = p.Results.debug;
tol = p.Results.tol;
ttol = p.Results.ttol;
shift = p.Results.shift;
iterative_mult = p.Results.iterative_mult;
use_sinc = p.Results.use_sinc;
interval_report = p.Results.interval_report;
conditional_indices = p.Results.conditional_indices;
x0 = p.Results.x0;
anderson = p.Results.anderson;

switch fun
	case 'inv'
		[m, time] = eval_inv(pi0, rewards, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, anderson);
	case 'inv2'
		[~, time1, y] = eval_inv(pi0, rewards, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, anderson);
		[m, time2] = eval_inv(pi0, -y, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, anderson);
		time = time1 + time2;
    
    case 'tta_variance'
        [m1, time1, y] = eval_inv(pi0, rewards, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, anderson);
		[m2, time2] = eval_inv(pi0, -y, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, anderson);
        m = m2 - m1^2;
		time = time1 + time2;
	
	case 'cond_etta'
		% Conditional expected time to absorption: we select as
		% absorbing_state to condition the first in the matrix
		% absorbing_states
		[m, time] = eval_cond_etta(pi0, R, W, absorbing_states, ...
								   conditional_indices, ...
								   algorithm, debug, tol, ttol, shift, ...
								   iterative_mult, use_sinc, interval_report);
							   
		if debug
			fprintf('EVAL_MEASURE :: cond_etta :: measure = %e\n', m);
			fprintf('EVAL_MEASURE :: cond_etta :: time    = %f\n', time);
		end
		
	otherwise
		error('Unsupported measure');
end

end

