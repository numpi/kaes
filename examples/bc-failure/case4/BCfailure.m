function [R, W, M] = BCfailure(n, topology, eta, zeta, lambdaR, lambdaD, lambdaE, lambdaEP, l)
%BCFAILURE produces matrices R and W and markings M
%   INPUTS: n:        number of system components,
%           topology: adjacency matrix of the failure propagation graph
%           lambdaR:  ...
%           lambdaD:  component failure detection rate
%           lambdaE:  local error
%           lambdaEP: failure propagation rate
%
% The model constructed is described in the paper "New Method for 
% the Evaluation of (Conditional) Time To Failure and its Moments
% Trough Implicit Reward Structure", by G. Masetti, L. Robol, S. Chiaradonna,
% F. Di Giandomenico, submitted.

% Check dimension
if size(topology) ~= [n, n]
    fprintf('ERROR: topology is required to be a %dx%d matrix\n',n);
    return;
end

R = cell(1, n);
W = cell(2*n, n);
M = cell(1,n);

% ordering of places: W, E, D, C, B, useful for M{i}

for i = 1 : n
    if i <= l
        % TODO without TB
        R{i} = zeros(7,7);
        R{i}(1,2) = 2*lambdaE(i);
        R{i}(2,3) = lambdaE(i);
        R{i}(2,5) = eta*lambdaD(i);
        R{i}(3,6) = 2*eta*lambdaD(i);
        R{i}(5,1) = zeta*lambdaR(i);
        R{i}(5,8) = (1-zeta)*lambdaR(i);
        R{i}(5,6) = lambdaE(i);
        R{i}(6,2) = zeta*lambdaR(i);
        R{i}(6,7) = eta*lambdaD(i);
        R{i}(7,5) = zeta*lambdaR(i);
        M{i} = [2 0 0 0 0
                1 1 0 0 0
                0 2 0 0 0
                0 0 0 0 1  % absorbing
                1 0 1 0 0 
                0 1 1 0 0
                0 0 2 0 0];
    else        
        R{i} = zeros(12,12);
        R{i}(1,2) = 2*lambdaE(i);
        R{i}(2,3) = lambdaE(i);
        R{i}(2,5) = eta*lambdaD(i);
        R{i}(3,6) = 2*eta*lambdaD(i);
        R{i}(5,1) = zeta*lambdaR(i);
        R{i}(5,8) = (1-zeta)*lambdaR(i);
        R{i}(5,6) = lambdaE(i);
        R{i}(6,2) = zeta*lambdaR(i);
        R{i}(6,7) = eta*lambdaD(i);
        R{i}(7,5) = zeta*lambdaR(i);
        % additional edges wrt the other case
        R{i}(6,9) = (1-zeta)*lambdaR(i);
        R{i}(7,10) = 2*(1-zeta)*lambdaR(i);
        R{i}(8,9) = lambdaE(i);
        R{i}(9,10) = eta*lambdaD(i);
        R{i}(10,8) = zeta*lambdaR(i);
        R{i}(10,11) = (1-zeta)*lambdaR(i);
        M{i} = [2 0 0 0 0
                1 1 0 0 0
                0 2 0 0 0
                0 0 0 0 1  % absorbing
                1 0 1 0 0 
                0 1 1 0 0
                0 0 2 0 0
                1 0 0 1 0
                0 1 0 1 0
                0 0 1 1 0
                0 0 0 2 0  % absorbing
                0 0 0 1 1]; % absorbing
    end
end

% manage TEP
for i = 1 : n
    for j = 1 : n
        if topology(i,j)
            if j <= l
                W{i,j} = zeros(7,7);
            else
                W{i,j} = zeros(12,12);
            end
            
            if i == j
                W{i,j}(2,2) = lambdaEP(i);
                W{i,j}(3,3) = 2*lambdaEP(i);
                W{i,j}(6,6) = lambdaEP(i);
                if j>l
                    W{i,j}(9,9) = lambdaEP(i);
                end
            else
                if j <= l
                    W{i,j} = zeros(7,7);
                    W{i,j}(1,3) = 1;
                    W{i,j}(2,3) = 1;
                    W{i,j}(5,3) = 1;
                    W{i,j}(6,3) = 1;
                    W{i,j}(7,3) = 1;
                else
                    W{i,j} = zeros(12,12);
                    W{i,j}(1,3) = 1;
                    W{i,j}(2,3) = 1;
                    W{i,j}(5,3) = 1;
                    W{i,j}(6,3) = 1;
                    W{i,j}(7,3) = 1;
                    % additional edges wrt the other case
                    W{i,j}(8,9) = 1;
                    W{i,j}(10,9) = 1;
                end
            end
        else
            % if component j is subject to NO exterior attack
            if j <= l
                W{i,j} = eye(7,7);
            else
                W{i,j} = eye(12,12);
            end
        end
    end
end

% manage TC
for i = n+1 : 2*n
    for j = 1 : n
        if j<=l
            W{i,j} = zeros(7,7);
        else
            W{i,j} = zeros(12,12);
        end

        if i-n == j
            W{i,j}(1,4) = 2*(1-eta)*lambdaD(i-n);
            W{i,j}(2,4) = 2*(1-eta)*lambdaD(i-n);
            W{i,j}(3,4) = 2*(1-eta)*lambdaD(i-n);
            W{i,j}(5,4) = 2*(1-eta)*lambdaD(i-n);
            W{i,j}(6,4) = 2*(1-eta)*lambdaD(i-n);
            W{i,j}(7,4) = 2*(1-eta)*lambdaD(i-n);
            if j>l
                W{i,j}(8,12) = (1-eta)*lambdaD(i-n);
                W{i,j}(9,12) = (1-eta)*lambdaD(i-n);
                W{i,j}(10,12) = (1-eta)*lambdaD(i-n);
            end
        else
            W{i,j}(1,4) = 1;
            W{i,j}(2,4) = 1;
            W{i,j}(3,4) = 1;
            W{i,j}(5,4) = 1;
            W{i,j}(6,4) = 1;
            W{i,j}(7,4) = 1;
			if j>l
                W{i,j}(8,12) = 1;
                W{i,j}(9,12) = 1;
                W{i,j}(10,12) = 1;
            end
        end
    end
end

% convert Rs and Ws in TT-format
for i = 1 : n    
    R{i} = tt_matrix(R{i});
end

for i = 1 : 2*n
    for j = 1 : n
        W{i,j} = tt_matrix(W{i,j});
    end
end

end
