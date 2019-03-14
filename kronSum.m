function C = kronSum(varargin)
%KRONSUM computes the Kronecher sum C = kron(A,IB) + kron(IA,B)
%   IA is the indentity matrix of the same size of A,
%   similarly for IB.

if length(varargin) > 3
    C = kronSum( kronSum(varargin{1}, varargin{2}, varargin{end}), varargin{3:end} );
else
    A = varargin{1};
    B = varargin{2};
    tol = varargin{3};
    
    if isa(B, 'tt_matrix')
        IA = tt_eye(size(A, 1), length(size(A)) / 2);
        IB = tt_eye(size(B, 1), length(size(B)) / 2);
    else
        IA = speye(size(A));
        IB = speye(size(B));
    end

    if isa(B, 'tt_matrix')
        C = tkron(A,IB) + tkron(IA,B);
    else
        C = kron(A,IB) + kron(IA,B);
    end
    
    if isa(C, 'tt_matrix')
        C = round(C, tol);
    end
end

end

