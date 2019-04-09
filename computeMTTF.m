function [m, time] = computeMTTF(R, W, ttol, tol, method, debug)
%COMPUTEMTTF Compute the MTTF for the system genrated by Q = R + W.
%
% M = COMPUTEMTTF(R, W, TTOL, TOL, METHOD, DEBUG) computes the MTTF (Mean
%     Time To Failure) for the continuous time Markov chain with 
%     infinitesimal generator Q := R + W, where R, W are expected to be 
%     given in TT (tensor train) format as in the TT-Toolbox.
%
%     TTOL and TOL are the tolerances for the truncation in the TT format, 
%     and for the stopping of the iteration in the method of choice. 
%
%     METHOD selects the computational engine. Available methods are:
%      - 'dense' constructs the matrix Q explicitly, and computes the MTTF 
%            using the formula:
%               MTTF = - pi0' * Q(1:end,1:end-1) \ Q(1:end-1,end)
%
%         Note that this method is given for reference only, and does not
%         scale to large systems. 
%
%      - 'expsums' computes the MTTF solving a linear system with Q - S,
%         where S is a proper rank 1 correction, using a Neumann expansion
%         and exploiting exponential sums in the inner iteration. The
%         computations are performed in dense arithmetic. 
%
%      - 'ttexpsums' is mathematically the same algorithm as 'expsums', but
%         is performed in the TT format. 
%
%      - 'ttexpsumst' is similar to 'ttexpsums', but works with Q' instead
%         of Q. This has the advantage that forces the zero entries
%         corresponding to non-reachable states in the vector computed in
%         the Neumann iteration, keeping  the TT-rank lower.
%         
%      - 'ttexpsums2' is a quadratically convergent matrix iteration that,
%         after k steps, obtains the same result of 2^k steps of the Neumann
%         iteration. Since it stores matrices in TT-form, instead of
%         vectors, this approach requires more storage with respect to
%         'ttexpsums' or 'ttexpsumst', but can be considerably faster.
%
%      - 'spantree' determines all the reachable states by computing a
%         spanning tree of the graph induced by the Markov chain, and then
%         computes the MTTF using the same approach as 'dense', but on the
%         matrix Q restricted to the reachable states.
%
%    DEBUG, if set to true, enables the debug output. 

time = 0;

if ~exist('debug', 'var')
    debug = false;
end

maxsteps = inf;

% This value might be tweaked depending on the spectrum location of the
% matrix R. At the moment this is done using an heuristic after calling
% minmaxeig -- but this approach could be refined.
expn  = 4;

k = length(R);
n = arrayfun(@(i) size(R{i}, 1), 1 : k);

if ~exist('method', 'var')
    method = 'ttexpsums';
end

DeltapC = diagblocks(R, W);
Deltap = -kronSum(DeltapC{:}, ttol);

for i = 1 : size(W, 1)
    WsyncInside = W{i,1};
    for j = 2 : k
        WsyncInside = tkron(WsyncInside, W{i,j});
    end
    
    if i > 1
        Wsync = Wsync + WsyncInside;
    else
        Wsync = WsyncInside;
    end
    
    Wsync = round(Wsync, ttol);
end

if exist('Wsync', 'var')
    QQ = kronSum(R{:}, ttol) + Wsync;
else
    QQ = kronSum(R{:}, ttol);
end

Delta = -diag(QQ * tt_ones(n));

A1 = kronSum(R{:}, ttol);
if exist('Wsync', 'var')
    A2 = round(Wsync + (Delta - Deltap), ttol);
else
    A2 = round(Delta - Deltap, ttol);
end

D = Deltap;

pi0 = createpi0(n);
en = createen(n);
S = round(tt_matrix(kron(en, (A1 + A2) * en)), ttol);

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

fprintf('> dimension = %d, cond = %e \n', prod(n), cnd);

if cnd < 2
    expn = 4;
elseif cnd < 4
    expn = 8;
end

% MTTF Q
switch method
    case 'dense'
        tic;
        FD = full(D); FS = full(S); FA1 = full(A1); FA2 = full(A2);
        m = (FD + FA1) \ ones(prod(n), 1);
        y = m;
        j = 1;
        while j < maxsteps
            y = -(FD + FA1) \ ((FA2 - FS) * y);
            m = m + y;
            j = j + 1;
        end
        m = m / scl;
        t = toc;
        time = t;
        fprintf('m = %e (full), time = %f sec\n', -m(1), t);
            
    case 'expsums'

        % MTTF Q expsums
        tic;
        FD = full(D); FS = full(S); FA1 = full(A1); FA2 = full(A2);
        II = expsuminv(-FD - FA1, expn);
        m = -II * ones(n^k, 1);
        y = m;
        j = 1;
        while j < maxsteps
            y = II * ((FA2 - FS) * y);
            m = m + y;
            j = j + 1;
        end
        m = m / scl;
        t = toc;
        time = t;
        fprintf('m = %e (exp sums full), time = %f sec\n', -m(1), t);
        
    case 'ttexpsums'
        % MTTD Q expsums tt
        tic;
        DA = R;
        for j = 1 : k
            DA{j} = (- DA{j} + DeltapC{j}) / scl;
        end
        y0 = tt_ones(n,k) - en;
        m = ttexpsummldivide(DA, y0, expn, ttol);
        y = m;
        nrmY = norm(y);
        j = 1;
        while j < maxsteps
            y = ttexpsummldivide(DA, (A2 - S) * y, expn, ttol);
            m = round(m + y, ttol);
            oldnrmY = nrmY;
            nrmY = norm(y); nrmM = norm(m);
            rho = nrmY / oldnrmY;
            err = nrmY / (1 - rho) / norm(m);
            if debug
                fprintf('Step %d, Neumann residue ~ %e, norm(m) = %e, erank = %f, spectral radius ~ %e\n', ...
                j, nrmY, nrmM, erank(m), rho);
                fprintf('Measure estimate: %e (err. estimate = %e, est. upper bound = %e)\n', ...
                    dot(m, pi0) / scl, err, dot(m, pi0) * (1 + err) / scl);
            end

            if rho < 1 && err < tol
                break
            end
            
            j = j + 1;
        end
        m = m / scl;
        t = toc;
        time = t;
        fprintf('m = %e (exp sums tt), time = %f sec\n', dot(pi0, m), t);
        
 case 'ttexpsumst'
        % MTTD Q expsums tt
        tic;
        DA = R;
        for j = 1 : k
            DA{j} = (- DA{j} + DeltapC{j})' / scl;
        end
        r = tt_ones(n,k) - en;
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
            err = nrmY / (1 - rho) / norm(m);
            if debug
            fprintf('Step %d, Neumann residue ~ %e, norm(m) = %e, erank = %f, erank y = %f, spectral radius ~ %e\n', ...
                j, nrmY, nrmM, erank(m), erank(y), rho);
            fprintf('Measure estimate: %e (err. estimate = %e, est. upper bound = %e)\n', ...
                dot(m, r) / scl, err, dot(m, r) * (1 + err) / scl);
            end

            if rho < 1 && rho <= oldrho * (1 + 1e-8) && err < tol
                break
            end
            
            j = j + 1;
        end
        m = m / scl;
        t = toc;
        time = t;
        fprintf('m = %e (exp sums tt), time = %f sec\n', dot(r, m), t);        

        
    case 'ttexpsums2'
        % MTTD Q expsums tt
        tic;
        DA = R;
        for j = 1 : k
            DA{j} = (- DA{j} + DeltapC{j}) / scl;
        end

        y = -ttexpsummldivide(DA, tt_ones(n), expn, ttol);
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
        t = toc;
        fprintf('m = %e (exp sums tt), time = %f sec\n', m, t);
        time = t;
        
    case 'spantree'
        tic;
        [Q, ~, ~, ~] = infgen(R, W, cell(1,1), 1e-8, true);
        
        if numel(Q) > 5e13 && false
            error('Q is too large: n = %d', numel(Q));
        end

        G = digraph(abs(Q) > 1e-3, 'omitselfloops');
        G = shortestpathtree(G, 1);
        G = graph(G.Edges, G.Nodes);
        t = minspantree(G, 'Root', 1);
        idx = unique(t.Edges.EndNodes);
        fprintf('Number of reachable states: %d\n', length(idx));
        QQ = Q(idx, idx);
        m = -QQ(1:end-1,1:end-1) \ ones(length(idx)-1, 1); m = m(1);
        t = toc;
        time = t;
        fprintf('m = %e (span tree), time = %f sec\n', m, t);
        
    case 'gmres'
        % Construct R - 
        DeltapC = diagblocks(R, W, 0);
        DA = R;
        for j = 1 : k
            DA{j} = round(R{j} - DeltapC{j}, ttol);
        end
        
        QQ = kronSum(R{:}, ttol) + Wsync;
        Delta = diag(QQ * tt_ones(n));
        
        QQ = round(QQ - Delta, ttol);
        
        kk = prod(size(pi0));
        b = round(tt_ones(n) - en, ttol);
        %l = gmres(full(QQ), full(b), kk, ...
        %        tol, kk, full(kronSum(DA{:}, ttol)));
        %m = -dot(full(pi0), l);
        
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
        expinv = @(x) -ttexpsummldivide(DA, x, 12, ttol);
        
        tic;
        l = tt_gmres_block(...
            @(x,ttol) { round(expinv(QQ*x{1}), ttol) }, ...
            expinv(b), ttol, ...
            'restart', 100, 'max_iters', 100, 'tol_exit', tol);
        m = -dot(pi0, l{1});
        t = toc;
        
        if debug
            fprintf('m = %e (gmres), time = %f sec\n', m, t);
        end
        
        time = t;
        
        
    otherwise
        error('Unsupported method');
        
end

end

