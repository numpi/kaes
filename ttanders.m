function [y, DX, DF] = ttanders(x, gx, nv, beta, tol, DX, DF)
%TTANDERS Tensor Train Anderson acceleration

if ~exist('DX', 'var') || isempty(DX)
    DX = {};
    DF = {};
end

% Number of vectors that we have already
m = length(DX);

f = round(gx - x, tol);

if m == 0
    y = round(x + beta * f, tol);
else
    % Store in the last column the updated values of x, f
    DX{m} = round(x - DX{m}, tol);
    DF{m} = round(f - DF{m}, tol);

    % Solve the linear system gam = DF \ f, using regularized least squares;
    % unfortunately, we do not have an implemented QR least square solver in
    % the toolbox, so we do it from scratch
    R = zeros(m);
    ff = zeros(m, 1);

    A = DF;

    for j = 1 : m
        % Make the vector A(:,j) orthogonal to the others
        for i = 1 : j - 1
            R(i,j) = dot(A{i}, A{j});
            A{j} = round(A{j} - R(i,j) * A{i}, tol);
        end

        % Compute the norm of A(:,j)
        R(j,j) = norm(A{j});
        A{j} = A{j} / R(j,j);

        % Compute f(j)
        ff(j) = dot(A{j}, f);
    end

    % Regularized inverse
    gam = pinv(R, 1e-8) * ff;

    y = x;
    for j = 1 : m
        y = round(y - gam(j) * DX{j}, tol);
    end

    if (beta ~= 0.0)
        y = round(y + beta * f, tol);
        for j = 1 : m
            y = round(y - beta * DF{j} * gam(j), tol);
        end
    end
end

% Store the old vectors x, f for the next run
DX{end+1} = x;
DF{end+1} = f;

if m >= nv
    DX = DX(2:end);
    DF = DF(2:end);
end

end

