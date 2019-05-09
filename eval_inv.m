function [m, time] = eval_inv(pi0, r, R, W, absorbing_states, ...
							  algorithm, debug, tol, ttol)
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

% This value might be tweaked depending on the spectrum location of the
% matrix R. At the moment this is done using an heuristic after calling
% minmaxeig -- but this approach could be refined.
expn  = 4;

DeltapC = diagblocks(R, W);
Deltap = - round(ktt_kronsum(DeltapC{:}), ttol);

Wsync = ktt_zeros(n, n, fmt);
for i = 1 : size(W, 1)
	Wsync = round(Wsync + ktt_kron(W{i,:}), ttol);
end

QQ =round(ktt_kronsum(R{:}) + Wsync, ttol);

Delta = -diag(QQ * ktt_ones(n, fmt));

A1 = round(ktt_kronsum(R{:}), ttol);
A2 = round(Wsync + (Delta - Deltap), ttol);

D = Deltap;

% Construct the low-rank correction to implicitly deflated the absorbing states
S = ktt_zeros(n, n, fmt);
for i = 1 : size(absorbing_states, 1)
	ej = ktt_ej(n, absorbing_states(i,:), fmt);
	S = round(S + ktt_outerprod(A1*ej + A2*ej, ej), ttol);
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

if cnd < 2
    expn = 4;
elseif cnd < 4
    expn = 8;
end



switch algorithm
	
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
            X = round( X * X, ttol * nrmX0 / nrmX, maxrank );

            nrmX = norm(X);
			nrmX0 = max(nrmX, nrmX0);
            
            m = -dot(y, pi0) / scl;

            fprintf('Step %d, Residue ~ %e, erank = %f %f, m = %f\n', ...
                j, nrmX, erank(X), erank(y), -dot(y, pi0) / scl);

            if nrmX < 1e-8 || (m - oldm) < m * tol
                break;
            end
            
            j = j + 1;
        end

        m = - dot(y, pi0);
        m = m / scl;
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
		abs_idx = [];
		for j = 1 : size(absorbing_states, 1)
			nn = cumprod(n, 'reverse');
			abs_idx = [abs_idx, absorbing_states(j,end) + ...
				       sum((absorbing_states(j,1:end-1)-1).*nn(2:end))];
		end
		
		idx = setdiff(idx, abs_idx);
		m = -Q(idx,idx) \ fr(idx); m = dot(m, fpi0(idx));
		
        t = toc(tspantree);
        time = t;
        fprintf('m = %e (span tree), time = %f sec\n', m, t);
		
	otherwise
		error('Unsupported algorithm');
		
end



end

