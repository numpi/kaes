function [m, time, y] = eval_inv(pi0, r, R, W, absorbing_states, ...
							  algorithm, debug, tol, ttol, shift, ...
							  iterative_mult, use_sinc, ...
                              interval_report, x0, anderson)
%EVAL_INV 

k = length(R);
n = arrayfun(@(i) size(R{i}, 1), 1 : k);

if isempty(absorbing_states)
	warning('No absorbing states specified, assuming it is the last');
	% Unless there is a specific choice, we assume that there is a unique 
	% absorbing state, and it is the last one.
	absorbing_states = n;
end

fmt = ktt_format(R{1});

% Maximum number of steps for the iterative method of choice
maxsteps = inf;

DeltapC = diagblocks(R, W, shift);
Deltap = - round(ktt_kronsum(DeltapC{:}), ttol);

Wsync = ktt_zeros(n, n, fmt);
for i = 1 : size(W, 1)
	Wsync = round(Wsync + ktt_kron(W{i,:}), ttol);
end

QQ = round(ktt_kronsum(R{:}) + Wsync, ttol);

Delta = -round(diag(QQ * ktt_ones(n, fmt)), ttol);

A1 = round(ktt_kronsum(R{:}), ttol);
A2 = round(Wsync + (Delta - Deltap), ttol);

D = Deltap;

% Construct the low-rank correction to implicitly deflate the absorbing states
if isa(absorbing_states, 'tt_matrix')
	S = absorbing_states;
	S = round(A1 * S, ttol) + round(A2 * S, ttol);
	S = round(S, ttol);
else
	S = ktt_zeros(n, n, fmt);
	for i = 1 : size(absorbing_states, 1)
		ej = ktt_ej(n, absorbing_states(i,:), fmt);
		S = round(S + ktt_outerprod(A1*ej + A2*ej, ej), ttol);
	end
end

DA = R;
for j = 1 : k
    DA{j} = (- DA{j} + DeltapC{j});
end

% Compute minimum and maximum eigenvalues of a Kronecker sum
[scl, cnd] = minmaxeig(DA);

if scl ~= 0
    cnd = cnd / scl;
else
    cnd = 0;
    scl = 1;
end

D = D / scl;
A1 = A1 / scl;
A2 = A2 / scl;
S = S / scl;

if debug
	fprintf('> dimension = %d, cond = %e \n', prod(n), cnd);
end

if use_sinc
	expn = ceil(- 6 * log(ttol)/pi);
	
	if debug
		fprintf('> Selecting degree of exponential sums: %d\n', expn);
	end
else
	if cnd < 3
		expn = 4;
	elseif cnd < 4
		expn = 6;
	else
		expn = 8;
	end
end

if use_sinc
    expsums_method = 'sinc';
else
    expsums_method = 'lag';
end

switch algorithm
	
	case 'amen'
		timer = tic;
		
		Q = round(QQ + Delta - scl * S, ttol);
		xx = amen_block_solve({ Q }, { r }, ttol, 'nswp', 1000, 'tol_exit', tol);
		m = -dot(pi0, xx);
		
		res = norm(Q * xx - r) / norm(r);
		
		if res > tol
			error('AMEN did not converge within the prescribed number of sweeps');
		end
		
		time = toc(timer);
        y = xx;
		
		fprintf('m = %e (AMEn), time = %f sec\n', m, time);
        
    case 'ament'
        timer = tic;
        
        % We use R{:} as a preconditioner
        Qt = round(QQ + Delta - scl * S, ttol)';
        
        xx = amen_block_solve({ Qt }, { pi0 }, ttol, 'nswp', 1000, 'tol_exit', tol);
		m = -dot(r, xx);
		
		%res = norm(Qt * xx - rt) / norm(rt);
        
        %if res > tol
		%	error('AMEN did not converge within the prescribed number of sweeps');
        %end		
        
        time = toc(timer);
		fprintf('m = %e (AMEn + exp sums), time = %f sec\n', m, time);
        
    case 'tt-regular-splitting'
        tic;
        
        % DA are the factors of the Kronecker sum for (gamma * eye(n) - R)
        % Note that the sign is changed with respect to R - gamma*eye. 

        en = ktt_ej(n, n, 'tt');

        
        % In case we need to debug
        %M = full(ktt_kronsum(R{:})) + full(Deltap); M(:,end) = M(:,end) - full(ktt_kronsum(R{:})) * full(en);
        %S = zeros(prod(n)); S(:,end) = full(ktt_kronsum(R{:}) + Wsync + Delta) * full(en); S(end,end) = S(end,end) + 1;
        %N = -full(Delta) + full(Deltap) - full(Wsync) + S; N(:,end) = N(:,end) - full(ktt_kronsum(R{:})) * full(en);
        
        % Precompute f                
        % f = (R - gamma * eye(prod(n))) \ R(:,end);
        Rend = round(ktt_kronsum(R{:}) * en, ttol);
        f = ttexpsummldivide(DA, -Rend, expn, ttol, expsums_method);        
        g = f * (1./ (1 - dot(en, f)));
        
        Su = round(Rend + (Wsync*en) + (Delta*en), ttol);
                
        Mb = ttexpsummldivide(DA, -r, expn, ttol, expsums_method); 
        Mb = round(Mb - g * dot(en, Mb), ttol);
        
        if isempty(x0)
            x0 = Mb;
        end
        
        x = x0;
        
        nrmx0 = norm(Mb);
        
        j = 0;
        
        % Data for Anderson acceleration
        if anderson        
            % Q is only needed in case we want to compute the residual
            % during Anderson acceleration. 
            Q = infgen(R, W, tol, false);
            DXX = []; DFF = []; xa = x; nrmr = norm(r);
        end
        
        while j < maxsteps
            j = j + 1;
            
            % Adjust the truncation tolerance based on the accuracy that we
            % have achieved as of now
            ltol = ttol; % max(ttol, rel_change / cnd / 10);
            
            xold = x;
            % x = M \ (N * x) + x0;
            
            % Compute l = N * x
            l = -Wsync * x;
            l = round(l - Delta * x + Deltap * x, ltol);
            enx = dot(en, x);
            l = round(l + Su * enx - Rend * enx, ltol);
            
            % Solve the linear system M*x = l
            x = ttexpsummldivide(DA, -l, expn, ltol, expsums_method); 
            x = round(x - g * dot(en, x), ltol);
                        
            % Estimate the spectral radius            
            rho = norm(x) / norm(xold);
            
            % Update the iterate
            x = x + Mb;
            
            % Use Anderson acceleration, if enabled
            if anderson
                ka = 3;
                [xa, DXX, DFF] = ttanders(xold, x, ka, 1.0, ttol, DXX, DFF);            
                
                % Compute Anderson residual
                res_a = norm(Q * xa - r) / nrmr;
            else
                res_a = inf;
            end
            
            % Estimate for the error
            if rho < 1
                err_est = nrmx0 / norm(x) * rho^(j+1) / (1 - rho);
            else
                err_est = inf;
            end            
            
            if debug
                if anderson
                    fprintf('Step %d, err. est. = %e, m = %e (Anderson: %e / res = %e), ranks = %d, rho = %e \n', ...
                        j, err_est, -dot(pi0, x), -dot(pi0, xa), res_a, max(rank(x)), rho);
                else
                    fprintf('Step %d, err. est. = %e, m = %e, ranks = %d, rho = %e\n', ...
                        j, err_est, -dot(pi0, x), max(rank(x)), rho);
                end
            end
            
            if err_est < tol
                break;
            end
            
            if res_a < tol
                x = xa;
                break;
            end
        end
        
        m = -dot(pi0, x);
        
        y = x;
        
        t = toc;
        time = t;
        fprintf('m = %e (tt-regular-splitting), time = %f sec\n', m(1), t);       
        
  case 'tt-regular-splitting2'
        tic;
        
        % DA are the factors of the Kronecker sum for (gamma * eye(n) - R)
        % Note that the sign is changed with respect to R - gamma*eye. 

        en = ktt_ej(n, n, 'tt');
        
        % In case we need to debug
        %M = full(ktt_kronsum(R{:})) + full(Deltap); M(:,end) = M(:,end) - full(ktt_kronsum(R{:})) * full(en);
        %S = zeros(prod(n)); S(:,end) = full(ktt_kronsum(R{:}) + Wsync + Delta) * full(en); S(end,end) = S(end,end) + 1;
        %N = -full(Delta) + full(Deltap) - full(Wsync) + S; N(:,end) = N(:,end) - full(ktt_kronsum(R{:})) * full(en);
        
        % Precompute f                
        % f = (R - gamma * eye(prod(n))) \ R(:,end);
        Rend = round(ktt_kronsum(R{:}) * en, ttol);
        f = ttexpsummldivide(DA, -Rend, expn, ttol, expsums_method);        
        g = f * (1./ (1 - dot(en, f)));
        
        % S = Su * Sv', Sv = en
        Su = round(Rend + (Wsync*en) + (Delta*en), ttol);
                
        % x0 = M \ full(r);
        % System to solve A \ r
        % x = (R - gamma * eye(prod(n))) \ r; 
        % x = x - g * x(end);
        x0 = ttexpsummldivide(DA, -r, expn, ttol, expsums_method); 
        x0 = round(x0 - g * dot(en, x0), ttol);
        
        x = x0;        
        nrmx0 = norm(x0);
        
        % Compute the iteration matrix P = M \ N
        N = ktt_outerprod(Su, en);
        N = round(N - ktt_outerprod(Rend, en), ttol);
        N = round(N - Delta + Deltap, ttol);
        N = round(N - Wsync, ttol);
        
        P = ttexpsummldivide(DA, -N, expn, ttol, expsums_method); 
        P = round(P - ktt_outerprod(g, P'*en), ttol);
        
        j = 0;
        
        while j < maxsteps
            j = j + 1;
            
            xold = x;
            
            x = round(P * x, ttol);
            
            % Estimate the spectral radius            
            rho = ( norm(x) / norm(xold) )^(2^(-j+1));
            
            x = round(xold + x, ttol);
            P = P * P;
            P = round(P, ttol);
                                    
            % Estimate for the error
            if rho < 1
                err_est = nrmx0 / norm(x) * rho^(2^j+1) / (1 - rho);
            else
                err_est = inf;
            end            
            
            if debug
                fprintf('Step %d, err. est. = %e, m = %e, ranks = %d (x), %d (P), rho = %e\n', ...
                    j, err_est, -dot(pi0, x), max(rank(x)), max(rank(P)), rho);
            end
            
            if err_est < tol
                break;
            end
        end
        
        m = -dot(pi0, x);
        
        t = toc;
        time = t;
        fprintf('m = %e (tt-regular-splitting2), time = %f sec\n', m(1), t);            
        
  case 'tt-regular-splitting-hybrid'
        tic;
        
        % DA are the factors of the Kronecker sum for (gamma * eye(n) - R)
        % Note that the sign is changed with respect to R - gamma*eye. 

        en = ktt_ej(n, n, 'tt');
        
        % In case we need to debug
        %M = full(ktt_kronsum(R{:})) + full(Deltap); M(:,end) = M(:,end) - full(ktt_kronsum(R{:})) * full(en);
        %S = zeros(prod(n)); S(:,end) = full(ktt_kronsum(R{:}) + Wsync + Delta) * full(en); S(end,end) = S(end,end) + 1;
        %N = -full(Delta) + full(Deltap) - full(Wsync) + S; N(:,end) = N(:,end) - full(ktt_kronsum(R{:})) * full(en);
        
        % Precompute f                
        % f = (R - gamma * eye(prod(n))) \ R(:,end);
        Rend = round(ktt_kronsum(R{:}) * en, ttol);
        f = ttexpsummldivide(DA, -Rend, expn, ttol, expsums_method);        
        g = f * (1./ (1 - dot(en, f)));
        
        % S = Su * Sv', Sv = en
        Su = round(Rend + (Wsync*en) + (Delta*en), ttol);
                
        % x0 = M \ full(r);
        % System to solve A \ r
        % x = (R - gamma * eye(prod(n))) \ r; 
        % x = x - g * x(end);
        x0 = ttexpsummldivide(DA, -r, expn, ttol, expsums_method); 
        x0 = round(x0 - g * dot(en, x0), ttol);
        
        x = x0;        
        nrmx0 = norm(x0);
        
        % Compute the iteration matrix P = M \ N
        N = ktt_outerprod(Su, en);
        N = round(N - ktt_outerprod(Rend, en), ttol);
        N = round(N - Delta + Deltap, ttol);
        N = round(N - Wsync, ttol);
        
        P = ttexpsummldivide(DA, -N, expn, ttol, expsums_method); 
        P = round(P - ktt_outerprod(g, P'*en), ttol);
        
        j = 0;
        
        while j < maxsteps
            j = j + 1;
            
            xold = x;
            
            x = round(P * x, ttol);
            
            % Estimate the spectral radius            
            rho = ( norm(x) / norm(xold) )^(2^(-j+1));
            
            x = round(xold + x, ttol, 100);
            P = round(P * P, ttol, 10);
                                    
            % Estimate for the error
            if rho < 1
                err_est = nrmx0 / norm(x) * rho^(2^j+1) / (1 - rho);
            else
                err_est = inf;
            end            
            
            if debug
                fprintf('Step %d, err. est. = %e, m = %e, ranks = %d (x), %d (P), rho = %e\n', ...
                    j, err_est, -dot(pi0, x), max(rank(x)), max(rank(P)), rho);
            end
            
            if err_est < tol
                break;
            end
        end
        
        % m = -dot(pi0, x);
        
        % Refine with the standard iteration
        m = eval_inv(pi0, r, R, W, absorbing_states, ...
                  'tt-regular-splitting', debug, tol, ttol, shift, ...
                  iterative_mult, use_sinc, interval_report, x);
        
        t = toc;
        time = t;
        fprintf('m = %e (tt-regular-splitting-hybrid), time = %f sec\n', m(1), t);        
        
  case 'tt-regular-splitting-iterative-refinement'
        tic;
        
        % DA are the factors of the Kronecker sum for (gamma * eye(n) - R)
        % Note that the sign is changed with respect to R - gamma*eye. 

        en = ktt_ej(n, n, 'tt');
        
        % In case we need to debug
        %M = full(ktt_kronsum(R{:})) + full(Deltap); M(:,end) = M(:,end) - full(ktt_kronsum(R{:})) * full(en);
        %S = zeros(prod(n)); S(:,end) = full(ktt_kronsum(R{:}) + Wsync + Delta) * full(en); S(end,end) = S(end,end) + 1;
        %N = -full(Delta) + full(Deltap) - full(Wsync) + S; N(:,end) = N(:,end) - full(ktt_kronsum(R{:})) * full(en);
        
        % Precompute f                
        % f = (R - gamma * eye(prod(n))) \ R(:,end);
        Rend = round(ktt_kronsum(R{:}) * en, ttol);
        f = ttexpsummldivide(DA, -Rend, expn, ttol, expsums_method);        
        g = f * (1./ (1 - dot(en, f)));
        
        % S = Su * Sv', Sv = en
        Su = round(Rend + (Wsync*en) + (Delta*en), ttol);
        
        b = r;
        xf = ktt_zeros(n);
        
        %M = full(ktt_kronsum(R{:})) + full(Deltap); M(:,end) = M(:,end) - full(ktt_kronsum(R{:})) * full(en);        
        M = round(ktt_kronsum(R{:}) + Deltap - ktt_outerprod(ktt_kronsum(R{:}) * en, en), ttol);
        
        % Compute the iteration matrix P = M \ N
        N = ktt_outerprod(Su, en);
        N = round(N - ktt_outerprod(Rend, en), ttol);
        N = round(N - Delta + Deltap, ttol);
        N = round(N - Wsync, ttol);        
            
        P = ttexpsummldivide(DA, -N, expn, ttol, expsums_method); 
        P = round(P - ktt_outerprod(g, P'*en), ttol);        
        
        Pstart = P;
        
        ltol = ttol;
        ltol = 1e-4; % tol;
            
        for jj = 1 : 10
                
            % x0 = M \ full(r);
            % System to solve A \ r
            % x = (R - gamma * eye(prod(n))) \ r; 
            % x = x - g * x(end);
            x0 = ttexpsummldivide(DA, -b, expn, ltol, expsums_method); 
            x0 = round(x0 - g * dot(en, x0), ltol);
            x = x0;        
            nrmx0 = norm(x0);
            
            P = Pstart;
                       
            j = 0;

            while j < maxsteps
                j = j + 1;

                xold = x;

                x = round(P * x, ltol);

                % Estimate the spectral radius            
                nrmred = norm(x) / norm(xold);
                rho = ( nrmred )^(2^(-j+1));
                
                ltol = ttol; % min(tol, max(ttol, 1 - nrmred))

                x = round(xold + x, ltol);
                P = round(P * P, ltol);

                % Estimate for the error
                if rho < 1
                    err_est = nrmx0 / norm(x) * rho^(2^j+1) / (1 - rho);
                else
                    err_est = inf;
                end            

                if debug
                    fprintf('Step %d, err. est. = %e, m = %e, ranks = %d (x), %d (P), rho = %e\n', ...
                        j, err_est, -dot(pi0, x) - dot(pi0, xf), max(rank(x)), max(rank(P)), rho);
                end

                if err_est < tol
                    break;
                end
            end
            
            % Compute new b and update x
            if jj == 1
                xf = x;
            else
                xf = round(xf + x, ttol);
            end
            
            b = round(r - M * xf + N * xf + ktt_outerprod(Rend, en) * xf + ktt_outerprod(en, en) * xf, ltol);
            
            if norm(x) < norm(xf) * tol
                break;
            end
            
        end
        
        m = -dot(pi0, xf);
        
        t = toc;
        time = t;
        fprintf('m = %e (tt-regular-splitting-iterative-refinement), time = %f sec\n', m(1), t);        
        
    case 'dense-splitting'
        tic;
        
        % Q = full(infgen(R,W,ttol,true));
        Delta = full(Delta);
        R = full(ktt_kronsum(R{:}));
        W = full(Wsync);        
        Deltap = full(Deltap);        
        
        % gamma = max(1, norm(Delta, inf)) + 1;
        % Deltap = -gamma * eye(prod(n));
        
        S = zeros(prod(n)); S(:,end) = R(:,end) + W(:,end) + Delta(:,end);
        S(end,end) = S(end,end) + 1;
        
        % Start the iteration
        % M = R - gamma * eye(prod(n)); M(:,end) = M(:,end) - R(:,end);
        M = R + Deltap; M(:,end) = M(:,end) - R(:,end);
        % N = -Delta - gamma * eye(prod(n)) - W + S; N(:,end) = N(:,end) - R(:,end);
        N = -Delta + Deltap - W + S; N(:,end) = N(:,end) - R(:,end);
        
        fprintf('spectral radious=1-%e\n',1-max(abs(eig(M\N))));
        
        % Precompute f
        f = (R + Deltap) \ R(:,end);
        g = f ./ (1 - f(end));
        
        % S = Su * Sv', Sv = en
        Su = R(:,end) + W(:,end) + Delta(:,end);
        
        % System to solver A \ r
        x0 = M \ full(r);
        x = x0;
        
        j = 0;
        
        while j < maxsteps
            j = j + 1;
            
            xold = x;
            % x = M \ (N * x) + x0;
            
            % Compute l = N * x
            l = -W * x;
            % l = l - Delta * x - gamma * x;
            l = l - Delta * x + Deltap * x;
            l = l + Su * x(end);
            
            % Solve the linear system M*x = l
            % x = (R - gamma * eye(prod(n))) \ l; 
            x = (R + Deltap) \ l; 
            x = x - g * x(end);    
            
            % Update the iterate
            x = x + x0;
            
            if debug
                fprintf('Step %d, rel. change %e, m = %e\n', ...
                    j, norm(xold - x, 1), -dot(full(pi0), x));
            end
            
            if norm(xold - x, 1) < tol
                break;
            end
        end
        
        m = -dot(full(pi0), x);
        
        t = toc;
        time = t;
        fprintf('m = %e (dense-splitting), time = %f sec\n', m(1), t);        
		
	case 'dmrg'
		timer = tic;
		
		Q = round(QQ + Delta + scl * S, ttol);
		
		% Tolerance has been experimentally adjusted to deliver reasonable
		% results. 
		maxswp = 100;
		xx = dmrg_solve2(Q, r, tol * 1e-2, 'nswp', maxswp);
		
		res = norm(Q * xx - r) / norm(r);
		
		if res > tol
			error('DMRG did not converge within the prescribed number of sweeps');
		end
		
		m = -dot(pi0, xx);
		
		time = toc(timer);
		fprintf('m = %e (dmrg), time = %f sec\n', m, time);        
	
	 case 'ttexpsumst'
        % MTTD Q expsums tt
        timer = tic;
        DA = R;
        for j = 1 : k
            DA{j} = (- DA{j} + DeltapC{j})' / scl;
        end
        % r = ktt_ones(n) - en;
        y0 = pi0;
        m = ttexpsummldivide(DA, y0, expn, ttol);
        y = m;
        nrmY = norm(y);
        j = 1;
        rho = 1;
        while j < maxsteps
            y = ttexpsummldivide(DA, (A2 - S)' * y, expn, ttol);
            m = round(m + y, ttol);
            oldnrmY = nrmY;
            nrmY = norm(y); nrmM = norm(m);
            oldrho = rho;
            rho = nrmY / oldnrmY;    
            err = nrmY / (1 - rho) / nrmM;
            if debug && mod(j, interval_report) == 0
				fprintf('Step %d, Neumann residue ~ %e, norm(m) = %e, rank = %f, rank y = %f, spectral radius ~ %e\n', ...
					j, nrmY, nrmM, max(rank(m)), max(rank(y)), rho);
				fprintf('Measure estimate: %e (err. estimate = %e, est. upper bound = %e)\n', ...
					dot(m, r) / scl, err, dot(m, r) * (1 + err) / scl);
            end

            if rho < 1 && rho <= oldrho * (1 + 1e-2) && err < tol
                break
            end
            
            j = j + 1;
		end
		y = y / scl;
        m = m / scl;
        t = toc(timer);
        time = t;
        fprintf('m = %e (exp sums tt), time = %f sec\n', dot(r, m), t);   
		y = m;
		m = dot(r, m);
        
    case 'ttexpsums2'
        % MTTD Q expsums tt
        timer = tic;
        DA = R;
        for j = 1 : k
            DA{j} = (- DA{j} + DeltapC{j}) / scl;
        end

        y = -ttexpsummldivide(DA, r, expn, ttol);
        X = ttexpsummldivide(DA, (A2 - S), expn, ttol);

        maxrank = inf;

        nrmX0 = norm(X);
        nrmX = nrmX0;

        j = 1;
        while j < maxsteps
            oldm = -dot(y, pi0) / scl;
            
            z = round( X * y, ttol, maxrank );
            y = round( y + z, ttol, maxrank ); clear('z');
						
			if iterative_mult
				X = ktt_iterative_mult(X, min(1e-1, ttol * nrmX0 / nrmX), debug);
			else
			    X = round( X * X, min(1e-2, ttol * nrmX0 / nrmX), maxrank );
			end

            nrmX = norm(X);
			nrmX0 = max(nrmX, nrmX0);
            
            m = -dot(y, pi0) / scl;

            fprintf('Step %d, Residue ~ %e, rank = %d %d, m = %f\n', ...
                j, nrmX, max(rank(X)), max(rank(y)), -dot(y, pi0) / scl);

            if nrmX < sqrt(tol) || (m - oldm) < m * tol
                break;
            end
            
            j = j + 1;
		end

		y = y / scl;
        m = - dot(y, pi0);
        t = toc(timer);
        fprintf('m = %e (exp sums tt), time = %f sec\n', m, t);
        time = t;
        
    case 'spantree'
        tspantree = tic;
		
		% Construction of the matrix Q in sparse format
		SR = cell(1, k); SW = cell(size(W,1), k);
		for i = 1 : k
			SR{i} = sparse(full(R{i}));
			for j = 1 : size(W, 1)
				SW{j,i} = sparse(full(W{j,i}));
			end
		end
		Q = ktt_kronsum(SR{:}); 
		for i = 1 : size(SW, 1)
			Q = Q + ktt_kron(SW{i,:});
		end
		d = Q * ktt_ones(n, 'sparse'); Q = Q - spdiags(d, 0, prod(n), prod(n));
		
		% Find the reachable states by walking on the graph
		G = digraph(abs(Q) > 1e-3, 'omitselfloops');
        G = shortestpathtree(G, 1);
        G = graph(G.Edges, G.Nodes);
        t = minspantree(G, 'Root', 1);
        idx = unique(t.Edges.EndNodes);
		
		if debug
			fprintf(' - Number of reachable states: %d\n', length(idx));
		end
		
		% Construct dense representation of the rewards and of the pi0
		fr = full(r); fpi0 = full(pi0);
		
		% Construct the set of indices of the absorbing states
		if isa(absorbing_states, 'tt_matrix')
			S = full(absorbing_states);
			abs_idx = find(sum(S, 2) > 1e-8);
		else
			abs_idx = [];
			for j = 1 : size(absorbing_states, 1)
				nn = cumprod(n, 'reverse');
				abs_idx = [abs_idx, absorbing_states(j,end) + ...
						   sum((absorbing_states(j,1:end-1)-1).*nn(2:end))];
			end
		end
		
		idx = setdiff(idx, abs_idx);
		y = Q(idx,idx) \ fr(idx); m = -dot(y, fpi0(idx));
		
		% Convert to tt_tensor for compatibility
		yy = zeros(length(fr), 1); yy(idx) = y;
		y = tt_tensor(yy, size(r));
		
        t = toc(tspantree);
        time = t;
        fprintf('m = %e (span tree), time = %f sec\n', m, t);
		
    case 'gmres'
        % Construct R  
        DeltapC = diagblocks(R, W, shift);
        %DeltapC = cell(1, k);
        % for j = 1 : k
        %     DeltapC{j} = tt_matrix( diag(full(R{j}) * ones(n(j), 1)) );
        % end
        DA = R;
        for j = 1 : k
            DA{j} = round(R{j} - DeltapC{j}, ttol);
        end
        
        QQ = round(ktt_kronsum(R{:}) + Wsync, ttol);
        Delta = round(diag(QQ * ktt_ones(n)), ttol);
        
        QQ = round(QQ - Delta, ttol);
        
        b = r;
%		kk = prod(n);
%        l = gmres(full(QQ), full(b), kk, ...
%                tol, kk, full(kronSum(DA{:}, ttol)));
 %       m = -dot(full(pi0), l);
		% keyboard
        
        % Compute minimum and maximum eigenvalues of a Kronecker sum
        [scl, cnd] = minmaxeig(DA);

        if scl ~= 0
            cnd = cnd / scl;
        else
            cnd = 0;
            scl = 1;
        end
        
        QQ = QQ / scl;
        for j = 1 : length(DA)
            DA{j} = -DA{j} / scl;
        end
        b = b / scl;
        
        % expinv = @(x) tt_tensor(reshape(full(-kronSum(DA{:}, ttol)) \ full(x), n));
        expinv = @(x) -ttexpsummldivide(DA, x, 4, ttol);
        %DeltapB = round(Delta - kronSum(DeltapC{:}, ttol), ttol);
        %expinv = @(x) gmres_preconditioner(DA, DeltapB, x);
        
        % A good preconditioner ?
        % M = round(kronSum(R{:}, ttol) - Delta, ttol) / scl;
        % MM = round(M' * M, tol);
        % expinv = @(x) amen_solve2(MM, M' * x, 1e-2);
        
        timer = tic;
        l = tt_gmres_block(...
            @(x,ttol) { round(expinv(QQ*x{1}), ttol) }, ...
            expinv(b), 1e-6, ...
            'restart', 1800, 'max_iters', 500, 'tol_exit', tol);
        m = -dot(pi0, l{1});
        t = toc(timer);
        
        if debug
            fprintf('m = %e (gmres), time = %f sec\n', m, t);
        end
        
        time = t;		
		
	otherwise
		error('Unsupported algorithm');
		
end



end

