function Q = ktt_infgen(R, W)
%KTT_INFGEN 

k = length(R);
n = arrayfun(@(i) size(R{i}, 1), 1 : k);

fmt = ktt_format(R{1});

Q = ktt_kronsum(R{:}); 

for i = 1 : size(W, 1)
	Q = Q + ktt_kron(W{i,:});
end

d = Q * ktt_ones(n, fmt);

switch fmt
	case 'sparse'
		Q = Q - spdiags(d, 0, prod(n), prod(n));
	case 'tt'
		Q = Q - diag(d);
end

end

