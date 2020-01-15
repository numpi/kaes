function Y = ttexpsummldivide(A, B, n, tol, method, rmax)
%TTEXPSUMINV Compute A \ B using exponential sums. 

if ~exist('method', 'var')
    method = 'lag';
end

if ~exist('rmax', 'var')
    rmax = inf;
end

[a,b] = expsums(n, method);
k = length(A);

II = tt_matrix(expm(-b(1) * full(A{1})));
for i = 2 : k
    II = tkron(tt_matrix(expm(-b(1) * full(A{i}))), II);
end

Y = a(1) * II * B;

for j = 2 : n
    II = tt_matrix(expm(-b(j) * full(A{1})));
    for i = 2 : k
        II = tkron(tt_matrix(expm(-b(j) * full(A{i}))), II);
    end

    Y = round(Y + a(j) * II * B, tol, rmax);
end


end

