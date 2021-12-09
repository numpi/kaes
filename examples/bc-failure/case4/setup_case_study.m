% Topology
rng(2)

% Tolerances
ttol = 1e-10;
tol = 1e-4;

%  topology = full((eye(n,n)+sprand(n,n, 1 / n)) > 0)
if ~exist('topology', 'var')
    topology = createTopology(n, 0.2, 'starnoloops');
end

% Topology generate with Erdos plus a superdiagonal correction equal to all
% ones. 
load('topology', 'T');
topology = T(1:n, 1:n);
% xx = rand(n-1,1) > .5;
% % yy = rand(n-1,1) > .5;
% topology = eye(n) + diag(xx .* ones(n-1,1),1) + diag(ones(n-1,1),-1);
% topology = topology > 0;

disp(topology);

% Number of systems with only one absorbing state
l = 0;% when there is a catastrofic failure we cannot consider the presence of
% some componets with just one absorbing state

% Construct Rs and Ws
[R, W] = BCfailure(n, topology, eta, zeta, lambdaR, lambdaD, ...
	lambdaE, lambdaEP, l);

P = zeros(12,12); 
P(4,4) = 1;
P(11,11) = 1;
P(12,12) = 1;

w = cell(1, n);
for j = 1 : n
    if j <= l
        w{j} = tt_matrix(P(1:7,1:7));
    else
        w{j} = tt_matrix(P);
    end
end
absorbing_states = ktt_kron(w{:});

pi0 = ktt_ej([ 7*ones(1,l), 12*ones(1,n-l) ], ones(1,n));

v = [ 0 ; 0 ; 0 ; 1 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 1 ; 1 ];
w = cell(1, n); for j = 1 : n
    if j <= l
        w{j} = tt_tensor(v(1:7));
    else
        w{j} = tt_tensor(v); 
    end
end

r = round(ktt_ones([ 7*ones(1,l), 12*ones(1,n-l) ]) - ktt_kron(w{:}), 1e-8);