function T = createTopology(n, p)

T = eye(n, n);

T = T + T(randperm(n), :);

for i = 1: n
  if rand < p
    % Create a second connection on this row
    j = randi(n);
    T(i,j) = T(i,j) + 1;
  end
end

T = (T > 0);

end
