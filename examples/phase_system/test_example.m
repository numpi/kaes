% function test_example(n)
%TEST_EXAMPLE 

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

W = cell(nphases * n, n + 1);
for j = 1 : nphases
	S = zeros(nphases); S(j,j) = 1;
	for i = 1 : n
		W{i+(j-1)*n, 1} = S;
		for k = 1 : n
			W{i+(j-1)*n, k+1} = WW{j}{i,k};
		end
	end
end

Q = infgen(R, W, M, 1e-8, false);



% end

