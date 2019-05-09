% function test_example(n)
%TEST_EXAMPLE 

debug = true;
ttol = 1e-8;
tol  = 1e-4;

rng('default');

% Number of components
n = 2;

% Number of phases in the system
nphases = 3;

% Topology of failure propagation between the components
T = createTopology(n, 1/n, 'randperm');

% Topology and rates of the phases
TP = [ 0 1 0 ; ...
	   0 0 1 ; ...
	   1 0 0   ];
   
% We assume that the local transitions are phase independent
iso = rand(1, n);

WW = cell(1, nphases);

for j = 1 : nphases
	cas = rand(1, n);
	mu  = zeros(1, n);
	
	[R, WW{j}, M] = CCFModel(n, T, iso, cas, mu);
	% Q = infgen(R, W, M, 1e-8, false);
end

% We build the R and W factors of the large system, including the phases 
% as the first component. 
R = { TP, R{:} };

W = cell(1 + nphases * n, n + 1);
for j = 1 : nphases
	S = zeros(nphases); S(j,j) = 1;
	for i = 1 : n
		W{i+(j-1)*n, 1} = S;
		for k = 1 : n
			W{i+(j-1)*n, k+1} = WW{j}{i,k};
		end
	end
end

W{end,1} = -TP;
for j = 1 : n
	W{end,j+1} = [ 0 0 ; ...
		           0 1 ];
end

Q = infgen(R, W, M, 1e-8, false);

l = [1:3,5:7,9:11];
Qh = Q(l,l); 
w = - Qh \ ones(9,1);
fprintf('MTTF = %e\n', w(1));

% Convert everything to TT format
for i = 1 : length(R)
	R{i} = tt_matrix(R{i});
	for j = 1 : size(W, 1)
		W{j,i} = tt_matrix(W{j,i});
	end
end

Q = infgen(R, W, M, 1e-8, false);

absorbing_states = [ 1 : nphases ; 2 * ones(n, nphases) ];

m = invQh_reward(R, W, createpi0(2, n), tt_ones(2, n), absorbing_states, ttol, tol, debug);


% end

