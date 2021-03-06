function time = input_n_density01(n, method, topology)
%
% Constructs the example decribed in the paper "Stochastic modeling and 
%   evaluation of large interdependent composed models through Kronecker 
%   algebra and Exponential sums", by G. Masetti, L. Robol, S. Chiaradonna,
%   F. Di Giandomenico, submitted.

%  topology = full((eye(n,n)+sprand(n,n, 1 / n)) > 0)
if ~exist('topology', 'var')
    topology = createTopology(n, 0.2, 'starnoloops');
end

disp(topology);

  pp = symrcm(topology);

  lambda_min = 0.5;
  lambda_max = 1.5;
  lambda = lambda_min + (lambda_max-lambda_min)*rand(n,1)

  mu_min = 2000;
  mu_max = 3000;
  mu = mu_min + (mu_max-mu_min)*rand(n,1)

  p_min = 0.95;
  p_max = 1.0;
  p = p_min + (p_max-p_min)*rand(n,1)

  % Construct Rs and Ws
  [R, W] = largeModel(n, topology(pp,pp), lambda(pp), mu(pp), p(pp));

  % Compute the measure
  if n <= 9
    m = computeMTTF(R, W, 1e-6, 1e-3, 'spantree');
  end

  if strcmp(method, 'ttexpsums2')
    [m, time] = computeMTTF(R, W, 1e-8, 1e-3, method, true);
  else
    [m, time] = computeMTTF(R, W, 1e-8, 1e-3, method, true);
  end

end

