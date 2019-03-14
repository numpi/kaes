function [Q, M, Wsync, Rsys] = infgen(R, W, MM, tol, sp)
%INFGEN produces the infinitesimal generator matrix
%   Inputs: local rates matrices within cells R
%           synchronization matrices within cells W
%   Output: infinitesimal generator matrix
%
% If sp is true the output is a sparse matrix, otherwise is given in TT
% format. 

if sp
    Rloc = sparse(full(R{1,1}));
else
    Rloc = R{1,1};
end

% all element of R and W are assumed to be square matrices and 
% have the same size k
n = length(R);
k = arrayfun(@(i) size(R{i}, 1), 1:n);

M = MM{1};

for j = 2 : n
    if sp
        Rloc = kronSum(sparse(full(R{1,j})), Rloc, tol);
    else
        Rloc = kronSum(Rloc, R{1,j}, tol);
    end
end

for i = 1 : size(W, 1)
    if sp
        WsyncInside = sparse(full(W{i,1}));
    else
        WsyncInside = W{i,1};
    end
    
    for j = 2 : n
        if isa(WsyncInside, 'tt_matrix')
            WsyncInside = tkron(WsyncInside, W{i,j});    
        else
            WsyncInside = kron(sparse(full(W{i,j})), WsyncInside);
        end
    end
    
    if exist('Wsync', 'var')    
        Wsync = Wsync + WsyncInside;
    else
        Wsync = WsyncInside;
    end
    
    if isa(Wsync, 'tt_matrix')
        Wsync = round(Wsync, tol);
    end
end

if isa(Rloc, 'tt_matrix')
    if exist('Wsync', 'var')
        Q = round(Rloc + Wsync, tol);
    else
        Q = Rloc;
    end
else
    if exist('Wsync', 'var')
        Q = Rloc + Wsync;
    else
        Q = Rloc;
    end
end

if ~exist('Wsync', 'var')
    Wsync = 0;
end

if isa(Q, 'tt_matrix')
    e = tt_ones(k);
    s = Q * e;
    Q = Q - diag(s);
    Rsys = Rloc - diag(s);
else
    e = ones(size(Q, 2), 1);
    s = Q * e;
    Q = Q - spdiags(s, 0, size(Q, 1), size(Q, 2));
    Rsys = Rloc - spdiags(s, 0, size(Q, 1), size(Q, 2));
end

