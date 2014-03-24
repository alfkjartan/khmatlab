function Ainv = generalized_inverse(U,S,V) 
%   Ainv = generalized_inverse(U,S,V) 
% Returns the generalized inverse of A given its SVD.

% Kjartan Halvorsen
% 2006-02-15
  
tol= max(size(U)) * norm(S) * eps(class(S));

Sinv = zeros(size(S'));
for i=1:size(S,2)
  if S(i,i)> tol
    Sinv(i,i) = 1/S(i,i);
  end
end

Ainv = V*Sinv*U';