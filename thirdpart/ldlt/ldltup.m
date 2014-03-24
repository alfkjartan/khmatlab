%
% [newL,newD]=ldltup(L,D,v)
%
% This function updates the Cholesky factorization of A to the Cholesky 
% factorization of A+vv'.  i.e. If
%
%    A=L*D*L'
%
% then
%
%    A+v*v'=newL*newD*newL'
%
% It is assumed that A is symmetric and postive definite.
%
% Reference: Gill, Murray, and Wright, "Practical Optimization", p43.
% Author: Brian Borchers (borchers@nmt.edu)
%
function [newL,newD]=ldltup(L,D,v)
%
%  First, find the size of the matrix.  
%
n=size(L,1);
%
%  Initialize newL, newD, and t.
%
newL=L;
newD=D;
oldt=1;
%
%  The main loop.  See Gill, Murray, and Wright for details.  
%
for j=1:n,
  p=v(j);
  t=oldt+p^2/D(j,j);
  newD(j,j)=D(j,j)*t/oldt;
  beta=p/(D(j,j)*t);
  if (j < n),
    v(j+1:n)=v(j+1:n)-p*L(j+1:n,j);
    newL(j+1:n,j)=L(j+1:n,j)+beta*v(j+1:n);
  end;
  oldt=t;
end;
