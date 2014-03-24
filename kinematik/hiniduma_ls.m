function [m,aor]=hiniduma_ls(ys)
% [m,aor]=hiniduma_ls(ys)
% Least squares solution to CoR estimation according to Hiniduma et
% al. 
%
% Input
%    ys   ->   Observations. (3 x N x nmarkers)
% Output
%    m    <-   CoR estimation
%    aor  <-   AoR estimation

% Kjartan Halvorsen
% 2002-04-29
  
thr=1e-4;

% Compute the center of rotation
  
  [tre,N,nmarkers]=size(ys);
  
  A=zeros(3,3);
  b=zeros(3,1);
  
  
  for p=1:nmarkers
    v2=0;
    v3=zeros(3,1);
    vbar=(mean(ys(:,:,p)'))';
    vcov=cov(ys(:,:,p)',1)+vbar*vbar';
    
    for k=1:N
      v3=v3 + ys(:,k,p)*ys(:,k,p)'*ys(:,k,p);
      v2=v2 + ys(:,k,p)'*ys(:,k,p);
    end
    v3=v3/N;
    v2=v2/N;
  
    %    A=A+2*(vcov-vbar*vbar');
    A=A+2*cov(ys(:,:,p)',1);
    b=b+v3-vbar*v2;
    
  end

[U,S,V]=svd(A);

if S(3,3)<thr
  Ap=pinv(A);
  m=Ap*b;
else
  m=A\b;
end

aor=V(:,3);

