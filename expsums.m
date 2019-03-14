function [a, b, f] = expsums(n)
%EXPSUMS Approximate the inverse through exponential sums
%
% 1 / X ~ A(1) * EXP(-B(1) * X) + ... + A(N) * EXP(-B(N) * X), 
%
% where Re(X) > 0. 

[x, w] = lagpts(n);

a = w;
b = x;

for j = 1 : length(a)
    a(j) = a(j) * exp(b(j));
end

f = @(x) expsum_evaluate(a, b, x);

    function r = expsum_evaluate(a, b, x)
        r = a(1) * expm(-b(1) * x);
        for jj = 2 : length(b)
            r = r + a(jj) * expm(-b(jj) * x);
        end
    end

end

