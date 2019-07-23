function [m, time, y] = eval_inv(pi0, r, R, W, absorbing_states, ...
							  algorithm, debug, tol, ttol, shift, ...
							  iterative_mult, use_sinc)
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

switch algorithm
	
	case 'amen'
		timer = tic;
		
		Q = round(QQ + Delta + scl * S, ttol);
		xx = amen_block_solve({ Q.' * Q }, { Q.' * r }, tol^2, 'nswp', 1000);
		m = -dot(pi0, xx);
		
		res = norm(Q * xx - r) / norm(r);
		
		if res > tol
			error('AMEN did not converge within the prescribed number of sweeps');
		end
		
		time = toc(timer);
		
		fprintf('m = %e (AMEn), time = %f sec\n', m, time);
		
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
            if debug && mod(j, 50) == 0
				fprintf('Step %d, Neumann residue ~ %e, norm(m) = %e, erank = %f, erank y = %f, spectral radius ~ %e\n', ...
					j, nrmY, nrmM, erank(m), erank(y), rho);
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
        % Construct R - 
        DeltapC = diagblocks(R, W, 0);
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
		kk = prod(n);
        l = gmres(full(QQ), full(b), kk, ...
                tol, kk, full(kronSum(DA{:}, ttol)));
        m = -dot(full(pi0), l);
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
        expinv = @(x) -ttexpsummldivide(DA, x, 8, ttol);
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

