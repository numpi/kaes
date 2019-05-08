function [R, W, M] = CCFModel(n, topology, iso, cas, mu)
%CCFMODEL Common-Cause-of-Failure
%
% This model describes a system of N components, each can fail
% independently with rate ISO(j), j = 1, ..., N; alternatively, the
% component J can cause a failure to itself and all the components connected 
% according the given TOPOLOGY. 
%
% Each failed component is repaired with rate MU(j). 

R = cell(1, n);
W = cell(n, n);
M = cell(1,n);

for j = 1 : n
    R{j} = zeros(2,2);
    R{j}(1,2) = iso(j);
    R{j}(2,1) = mu(j);
    M{j} = [1 ; 0];
end

for i = 1 : n
   for j = 1 : n
      if ( topology(i,j) == 1 )
          W{i,j} = zeros(2,2);
          W{i,j}(1,2) = cas(i);
		  
		  % This synchronization transition is enabled only if the
		  % component J is not failed. 
		  if i ~= j
			W{i,j}(2,2) = 1.0;
		  end
      else
          % no interaction between i and j
    	  W{i,j} = eye(2,2);
      end
   end
end


end
