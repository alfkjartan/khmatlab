function [xhat,Phat,zhat,Pbar]=ekf_mechanism_fixlag(y,x0,P0,R1,R2,dt,tws,p0,m)
% [xhat,Phat,zhat,Pbar]=ekf_smooth_am(y,x0,P0,R1,R2,m)
%		         
% Extended fixed-lag smoothing Kalman Filter for tracking the state 
% of the system. 
%
% Input
%    y             ->   The data vector.
%    x0            ->   The initial state estimate
%    P0            ->   The initial state covariance
%			Name of an m-function 
%    R1            ->   The process noise covariance matrix
%    R2            ->   The measurement noise covariance matrix
%    m             ->   The lag.
% Output
%    xhat          <-   The state estimates
%    Phat          <-   The state covariance estimates
%

% Kjartan Halvorsen
% 2002-03-06

[d,N]=size(y);

n=length(x0);

n2=n/2;
Ft=[eye(n2) dt*eye(n2) ; zeros(n2,n2) eye(n2)];

if (size(R1,3)==1) % Process noise is time invariant
  isR1ti=1;
  % diagonalize
  if size(R1,2)==1
    G=diag(sqrt(R1));
  else    
    G=chol(R1)';
  end
else
  isR1ti=0;
end

%Gbar=kron([1;zeros(m,1)],zeros(n,size(R1,1)));
%Fbar=kron(diag(ones(m,1),-1),eye(n));

%Hbar=zeros(d,(m+1)*n);

zhat=zeros(n*(m+1),N);
Phat=zeros(n,n,N);

% so, the initialization
zt=[x0;zeros(n*m,1)];
Pt=zeros(n*(m+1),n*(m+1));
Pt(1:n,1:n)=P0;

for t=1:N
   % The observation function returns the estimates ouput yhat, and
   % the linearized observation function.
   [yhatt,Ht]=observe_mechanism_H(zt(1:n),y(:,t),tws,p0);

   %Hbar(1:d,1:n)=Ht;
   
   % The kalman gain
   if size(R2,3)>1   
      R2t=R2(:,:,t);
   else
      R2t=R2;
   end

%   size(Ht)
%   size(R2t)
%   n
   
   %HPHR=Hbar*Pt*Hbar'+R2t;
   HPHR=Ht*Pt(1:n,1:n)*Ht'+R2t;
  % if cond(HPHR)>1e12
      %Kt=Pt*Hbar'*psinv(HPHR);
  %    Kt=Pt(:,1:n)*Ht'*psinv(HPHR);
  % else
      %Kt=Pt*Hbar'*inv(HPHR);
      Kt=Pt(:,1:n)*Ht'*inv(HPHR);
  % end

   % The filter update
   ztt=zt + Kt*(y(:,t)-yhatt);
   zhat(:,t)=ztt;

   %Ptt=Pt - Kt*Hbar*Pt;
   Ptt=Pt - Kt*Ht*Pt(1:n,:);
   Phat(:,:,t)=Ptt(1:n,1:n);

   % Added 2002-05-14

   % Prediction
   xt=Ft*ztt(1:n);
   zt=[xt;ztt(1:m*n)];
   %Fbar(1:n,1:n)=Ft;

   % The process noise
   if ~isR1ti
     R1t=R1(:,:,t);
     % factorize
     G=chol(R1t)';
   end

   %Pt=Fbar*Ptt*Fbar'+Gbar*R1*Gbar';
   FtPtt=Ft*Ptt(1:n,1:end-n);
   Pt(1:n,1:n)=Ft*Ptt(1:n,1:n)*Ft';
   Pt(1:n,n+1:end)=FtPtt;
   Pt(n+1:end,1:n)=FtPtt';
   Pt(n+1:end,n+1:end)=Ptt(1:end-n,1:end-n);

   Pt(1:n,1:n)=Pt(1:n,1:n)+G*G';
end
	
Pbar=Phat;

xhat=[zhat(end-n+1:end,m+1:end) fliplr(reshape(ztt(1:m*n),n,m))];

% Reorganize the z matrix to a 3 dimensional array. The first layer
% is the filtered state estimates, the second the one step lag
% smoothed estimates, and so on so the m+1 layer is the lag m 
% smoothed estimates.

zorig=zhat;
zhat=zeros(n,N,m+1);

for l=1:m+1
   zhat(:,1:end-l+1,l)=zorig(n*(l-1)+1:n*l,l:end);
   if l>1
     zhat(:,end-l+2:end,l)=fliplr(reshape(ztt(1:(l-1)*n),n,l-1));
   end
end



