function [m,dir]=biascomp_ls(v)
% [m,dir]=biascomp_ls(v)
% Returns the bias-compensated least squares solution to the
% problem of estimating the center of rotation from the
% observations v.
%
% Input
%   v     ->  observations. (3 x N x nmarkers), where N is the
%             number of time samples, and nmarkers is the number of
%             markers.
% Output
%   m     <-  the estimated center of rotation (3 x 1) vector.
%   dir   <-  unit vector in the direction of the axis of rotation%

% Kjartan Halvorsen
% 2002-04-29
% 
% Revisions
% 2003-04-16  Added the direction of the axis of rotation as output.
% 2003-07-31  Marker data can contain zeros indicating missing
%             markers.
  
[tre,N,nmarkers]=size(v);

% Find indices of nonzero data
nonzero = cell(nmarkers,1);
for p=1:nmarkers
  nonzero{p} = find(~sum(v(:,:,p) == 0));
end

% Find A and b

b=zeros(3,1);
A=zeros(3,3);
  
for p=1:nmarkers
  vp = v(:,nonzero{p},p);
  M = size(vp,2);
  v2=0;
  v3=zeros(3,1);
  vbar = zeros(3,1); %(mean(v(:,:,p)'))';
  for k=1:M
      v3 = v3 + vp(:,k)*vp(:,k)'*vp(:,k);
      v2 = v2 + vp(:,k)'*vp(:,k);
  end
  
  vbar = mean(vp,2);
  v3 = v3/M;
  v2 = v2/M;
  
  A=A+2*cov(vp',1);
  b=b+v3-vbar*v2;
end

m=A\b;

dif=m;
thr=1e-14;

it=0;
while (norm(dif)>thr)
  it=it+1;
  
  % Find u, then ubar, and estimate of bias term
  u=zeros(size(v));
  ubar=zeros(3,nmarkers);
  for p=1:nmarkers
    u(:,:,p)=v(:,:,p)-kron(ones(1,N),m);
    ubar(:,p)=mean( u(:,nonzero{p},p), 2);
  end
  s2=sigma_estimate(u, nonzero);
  
  % Correct b, then run again
  db=zeros(3,1);
  for p=1:nmarkers
    M = length(nonzero{p});
    db=db+2*s2*(1-1/M)*ubar(:,p);
  end

  m_old=m;
  m=A\(b-db);
  
  dif=m-m_old;
  
end


% The axis of rotation
[V,D]=eig(A);
d=diag(D);
ind=find(d==min(d));
dir = V(:,ind(1));


function s2=sigma_estimate(u, nonzero)
%  s2=sigma_estimate(u, nonzero)
% Returns the estimate of the noise variance.
%
% Input
%    u     ->   the observations u=w+e. (3 x N x nmarkers)
%    nonzero->  cell array with indices of nonzero data
% Output
%    s2    <-   the noise variance estimate

% Kjartan Halvorsen
% 2002-04-29
  
[tre,N,nmarkers]=size(u);

% Find radius squared estimate u2.
u2=zeros(nmarkers,1);
for p=1:nmarkers
  up = u(:,nonzero{p},p);
  M = size(up,2);
  for k=1:M
    u2(p)=u2(p)+up(:,k)'*up(:,k);
  end
  u2(p)=u2(p)/M;
end

s2=0;
for p=1:nmarkers
  % The radius variance
  rvarp=0;
  up = u(:,nonzero{p},p);
  M = size(up,2);
  for k=1:M
    rvarp=rvarp+(up(:,k)'*up(:,k)-u2(p))^2;
  end
  s2=s2+1/(4*u2(p))*rvarp/M;
end

s2=s2/nmarkers;

