function Z = ktt_ones(varargin)
%KTT_ONES 

ninputs = length(varargin);

if ischar(varargin{end})
	fmt = varargin{end};
	ninputs = ninputs - 1;
else
	fmt = 'tt';
end

if ninputs > 2
	error('Unsupported number of inputs');
end

switch fmt
	case 'tt'
		if ninputs == 1
			Z = tt_ones(varargin{1}(end:-1:1));
		else
			if any(varargin{1} ~= varargin{2})
				error('ktt_ones only supports square TT matrices');
			end
			Z = tt_matrix(tt_ones(varargin{1}(end:-1:1) .* varargin{2}(end:-1:1)));
		end
	case 'sparse'
		m = varargin{1};
		if ninputs == 1
			Z = spones(prod(m), 1);
		else
			n = varargin{2};
			Z = spones(prod(m), prod(n));
		end
	otherwise
		error('Unsupported format');
end

end

