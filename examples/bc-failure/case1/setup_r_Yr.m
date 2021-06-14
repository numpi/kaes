% We need a different r for the _Yr cases. 
v = [ -2 ; 0 ; -1 ; 1 ; 1 ];
w = cell(1, n); 
for j = 1 : n
    if j <= l
        w{j} = tt_tensor(v(1:4)); 
    else
        w{j} = tt_tensor(v);
    end
end
r = round(ktt_ones([ 4*ones(1,l), 5*ones(1,n-l) ]) - ktt_kron(w{:}), 1e-8);