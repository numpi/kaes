function [x w v] = lagpts(n, int)
%LAGPTS Compute Laguerre integration points and weights by standard
%       Golub-Welsch algorithm. 
  
   alpha = 2*(1:n)-1;  beta = 1:n-1;     % 3-term recurrence coeffs
   T = diag(beta,1) + diag(alpha) + diag(beta,-1);  % Jacobi matrix
   [V,D] = eig(T);                       % eigenvalue decomposition
   [x,indx] = sort(diag(D));             % Laguerre points
   w = V(1,indx).^2;                     % Quadrature weights
   v = sqrt(x).*abs(V(1,indx)).';        % Barycentric weights
   v = v./max(v); v(2:2:n) = -v(2:2:n);
   
end
