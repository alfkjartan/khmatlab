function m=expm_rodrigues(what,theta)
% Calculates the exponential of what times theta using Rodrigues'
% formula.

% Kjartan Halvorsen
% 2001-09-28

tol=1e-12;

if nargin==1
  theta=1;
end

if (norm(what)==1)
   m = eye(3) + what*sin(theta) + what*what*(1-cos(theta));
elseif (norm(what)<tol)
   m=eye(3);
else
   nw=norm(what);
   m = eye(3) + what*sin(nw*theta)/nw + what*what*(1-cos(nw*theta))/nw^2;
end



