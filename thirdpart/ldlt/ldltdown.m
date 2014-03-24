%
% [newL,newD]=ldltdown(L,D,v)
%
% This function downdates the Cholesky factorization of A to the Cholesky 
% factorization of A-vv'.  i.e. If
%
%    A=L*D*L'
%
% then
%
%    A-v*v'=newL*newD*newL'
%
% It is assumed that A is symmetric and postive definite.  If A-v*v'
% would not be positive definite, then the diagonal elements will 
% be adjusted to make sure that D is > 0.  
%
% Reference: Gill, Murray, and Wright, "Practical Optimization", p43.
% Author: Brian Borchers (borchers@nmt.edu)
%
function [newL,newD]=ldltdown(L,D,v)
%
%  First, find the size of the matrix.  
%
n=size(D,1);
%
%  Initialize newL and newD.
%
newL=L;
newD=D;
%
%  Find the initial values of p and t.  
%
p=L\v;
oldt=1-p'*inv(D)*p;
%
%  Make sure that D is > 0.
%
if (oldt <= eps),
  oldt=eps;
end;
%
% The main loop.  See Gill, Murray, and Wright for details. 
%
for j=n:-1:1
  t=oldt+p(j)^2/D(j,j);
  newD(j,j)=D(j,j)*oldt/t;
  beta=-p(j)/(D(j,j)*oldt);
  v(j)=p(j);
  if (j < n),
    newL(j+1:n,j)=L(j+1:n,j)+beta*v(j+1:n);
    v(j+1:n)=v(j+1:n)+p(j)*L(j+1:n,j);
  end;
  oldt=t;
end;




