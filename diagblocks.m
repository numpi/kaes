function D = diagblocks(R, W, shift_param)
%DIAGBLOCKS Compute diagonal matrices D{j} such that 
%
%   kronSum(R{:}) <= kronSum(D{:})   component-wise

if ~exist('shift_param', 'var')
    shift_param = 1;
end

k = length(R);
n = arrayfun(@(i) size(R{i}, 1), 1 : k);

D = cell(1, k);

%for j = 1 : k
%    D{j} = tt_matrix(zeros(n(j)));
%end

nrm = 0;

for j = 1 : k
    D{j} = diag(R{j} * tt_ones(n(j), 1));
    nrm = nrm + max(abs(full(diag(D{j}))));
end

for j = 1 : k
    for i = 1 : size(W, 1)
        nrm = nrm + max(abs(full(diag(W{i,j}))));
    end
end 

for i = 1 : size(W, 1)
    % M = diag(W{i, 1} * tt_ones(n(1), 1));
    M = tt_eye(n(1)) * norm(full(W{i,1}), inf);
    for j = 2 : k
        M = M * norm(full(W{i,j}), inf);
    end
    %for j = 1 : k
    %    D{j} = D{j} + diag(W{i,j} * tt_ones(n(j), 1));
    %end
    % shift = 100 * max(norm(D{i}), norm(M)) * tt_matrix(eye(size(D{i}))) / size(W, 2);
    D{1} = round(D{1} + M, 1e-8);
end

for i = 1 : k
    shift = shift_param * nrm * tt_matrix(eye(size(D{i}))) / k;
    D{i} = D{i} + shift;
end


end

