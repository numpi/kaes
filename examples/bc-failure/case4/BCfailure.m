function [R, W, M] = BCfailure(n, topology, lambdaB, lambdaC, lambdaD, lambdaW, lambdaE, lambdaEP, l)
%BCFAILURE produces matrices R and W and markings M
%   INPUTS: n:        number of system components,
%           topology: adjacency matrix of the failure propagation graph
%           lambdaB:  benign component failure rates
%           lambdaC:  catastrofic component failure rate
%           lambdaD:  component failure detection rate
%           lambdaW:  component recovery rate
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

% ordering of places: W, E, D, C, B

for i = 1 : n
    if i <= l
        R{i} = zeros(7,7);
        R{i}(1,2) = 2*lambdaE(i);
        R{i}(2,3) = lambdaD(i);
        R{i}(2,4) = lambdaE(i);
        R{i}(4,5) = 2*lambdaD(i);
        R{i}(5,2) = lambdaW(i);
        R{i}(5,6) = lambdaD(i);
        R{i}(6,3) = 2*lambdaW(i);
        M{i} = [2 0 0 0 0
                1 1 0 0 0
                1 0 1 0 0
                0 2 0 0 0
                0 1 1 0 0
                0 0 2 0 0
                0 0 0 2 0];% absorbing
    else        
        R{i} = zeros(12,12);
        R{i}(1,2) = 2*lambdaE(i);
        R{i}(2,3) = lambdaD(i);
        R{i}(2,4) = lambdaE(i);
        R{i}(4,5) = 2*lambdaD(i);
        R{i}(5,2) = lambdaW(i);
        R{i}(5,6) = lambdaD(i);
        R{i}(6,3) = 2*lambdaW(i);
        % additional edges wrt the other case
        R{i}(3,8) = lambdaB(i);
        R{i}(5,9) = lambdaB(i);
        R{i}(8,9) = lambdaE(i);
        R{i}(9,10) = lambdaD(i);
        R{i}(10,8) = lambdaW(i);
        R{i}(10,12) = lambdaB(i);
        M{i} = [2 0 0 0 0
                1 1 0 0 0
                1 0 1 0 0
                0 2 0 0 0
                0 1 1 0 0
                0 0 2 0 0
                0 0 0 2 0 % absorbing
                1 0 0 0 1
                0 1 0 0 1
                0 0 1 0 1
                0 0 0 1 1 % absorbing
                0 0 0 0 2]; % absorbing
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
                W{i,j}(4,4) = lambdaEP(i);
                W{i,j}(5,5) = lambdaEP(i);
            else
                if j <= l
                    W{i,j} = zeros(7,7);
                    W{i,j}(1,2) = 1;
                    W{i,j}(2,4) = 1;
                    W{i,j}(3,4) = 1;
                    W{i,j}(5,4) = 1;
                    W{i,j}(6,4) = 1;
                else
                    W{i,j} = zeros(12,12);
                    W{i,j}(1,2) = 1;
                    W{i,j}(2,4) = 1;
                    W{i,j}(3,4) = 1;
                    W{i,j}(5,4) = 1;
                    W{i,j}(6,4) = 1;
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

        % TODO
        if i-n == j
            W{i,j}(1,7) = lambdaC(i-n);
            W{i,j}(2,7) = lambdaC(i-n);
            W{i,j}(3,7) = lambdaC(i-n);
            W{i,j}(4,7) = lambdaC(i-n);
            W{i,j}(5,7) = lambdaC(i-n);
            W{i,j}(6,7) = lambdaC(i-n);
            if j>l
                W{i,j}(8,11) = lambdaC(i-n);
                W{i,j}(9,11) = lambdaC(i-n);
                W{i,j}(10,11) = lambdaC(i-n);
            end
        else
            W{i,j}(1,7) = 1;
            W{i,j}(2,7) = 1;
            W{i,j}(3,7) = 1;
            W{i,j}(4,7) = 1;
            W{i,j}(5,7) = 1;
            W{i,j}(6,7) = 1;
			if j>l
                W{i,j}(8,11) = 1;
                W{i,j}(9,11) = 1;
                W{i,j}(10,11) = 1;
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
