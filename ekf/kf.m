function \
      [xhat,Phat]=kf(y,x0,P0,R1,R2,dt, H, \
                                                stepsahead)
%  [xhat,Phat]=kf(y,x0,P0,R1,R2,dt,H, stepsahead)
%		         
% k-steps-ahead Kalman filter. 
%
% Input
%    y             ->   The data vector.
%    x0            ->   The initial state estimate
%    P0            ->   The initial state covariance
%    R1            ->   The process noise covariance matrix
%    R2            ->   The measurement noise covariance matrix
%    dt            ->   The sampling time
%    H             ->   The measurement matrix (y = Hx)
%    stepsahead    ->   The number of time steps to look ahead
% Output
%    xhat          <-   The state estimates
%    Phat          <-   The state covariance estimates
%

% Kjartan Halvorsen
% 2011-12-02
%


[d,N]=size(y); % N is number of frames in data vector

n=length(x0); % Number of states in state vector

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
p=size(R1,1);

Gt=eye(n);

n2=n/2;
Ft=[eye(n2) dt*eye(n2) ; zeros(n2,n2) eye(n2)];
Ftinv=inv(Ft);

% Predictor matrix
Ftk = Ft^stepsahead;

if (size(R2,3)==1) % Measurement noise is time invariant
  isR2ti=1;
else
  isR2ti=0;
end

Ht = H; % Time invariant linear measurement function

% initialization

xt=x0;
Pt=P0;

xhat=zeros(n,N);
Phat=zeros(n,n,N);

restart = 0;
for t=1:N

%  warning('ekf-warning', ...
%          ['ekf_mechanism_predictor: t=', int2str(t)])


   % The measurement noise
   if isR2ti
     R2t=R2;
   else
     R2t=R2(:,:,t);
   end

   % The observations
   yhatt = Ht*xt;
   
   % The innovations
   et = y(:,t)-yhatt;

   % The kalman gain
   HPHR=Ht*Pt*Ht'+R2t;
   Kt=Pt*Ht'*inv(HPHR);

   % The filter update
   xtt=xt + Kt*et;

   %Ptt=Pt - Kt*Hbar*Pt;
   Ptt=Pt - Kt*Ht*Pt;
   
   % Prediction
   xt=Ft*xtt;

   % The process noise
   if ~isR1ti
     R1t=R1(:,:,t);
     % factorize
     G=chol(R1t)';
   end

   Pt=Ft*Ptt*Ft'+G*G';

  if (t==1)
    Ftkk = Ft;
    for (i=1:stepsahead)
      Ftkk = Ft*Ftkk;
      xhat(:,i) = Ftkk*xtt;
      Phat(:,:,i)= Ftkk*Pt*Ftkk'+G*G';
    end
  elseif(t<=N-stepsahead)
    xhat(:,t+stepsahead)=Ftk * xtt;
    Phat(:,:,t+stepsahead)= Ftk*Ptt*Ftk'+G*G';
  end


end
	
