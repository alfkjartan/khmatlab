function \
      [xhat,Phat]=ekf_mechanism_predictor(y,x0,P0,R1,R2,dt,tws,p0,rmse, \
                                                stepsahead)
%  [xhat,Phat]=ekf_mechanism_predictor(y,x0,P0,R1,R2,dt,tws,p0,rmse, stepsahead)
%		         
% k-steps-ahead Extended Kalman Filter. Modified from
% ekf_mechanism_fixint.
%
% Input
%    y             ->   The data vector.
%    x0            ->   The initial state estimate
%    P0            ->   The initial state covariance
%    R1            ->   The process noise covariance matrix
%    R2            ->   The measurement noise covariance matrix
%    dt            ->   The sampling time
%    tws           ->   Nested cell arrays of twists, determining
%                       the kinematic model
%    p0            ->   The reference positions of the markers
%    rmse          ->   The root mean square error in marker
%                       positions for the first frame. Used to
%                       detect divergence, and fall back to the
%                       inverse kinematics solution of function firststate.
%    stepsahead    ->   The number of time steps to look ahead
% Output
%    xhat          <-   The state estimates
%    Phat          <-   The state covariance estimates
%

% Kjartan Halvorsen
% 2011-09-12
%

rmse_tol = 1.3; % Tolerate 30% larger innovations than rmse of
                % initial marker residuals
rmse_tol = 2;  % Tolerate 200% larger innovations than rmse of
                % initial marker residuals

checkformarkererrors = 0; % Performs a check for errors in marker
                          % data at each time frame. 

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

   % The observation function returns the estimates ouput yhat, and
   % the linearized observation function.
   if checkformarkererrors
     [yhatt,Ht,slsk, y(:,t)]=observe_mechanism_H(xt,y(:,t),tws,p0,...
						 {Pt, R2t});
   else
     [yhatt,Ht]=observe_mechanism_H(xt,y(:,t),tws,p0);
   end
   
   % The innovations
   et = y(:,t)-yhatt;

   %etrms = et(find(y(:,t) ~= 0));
   %if ( sqrt(mean(etrms.*etrms)) > rmse_tol*rmse & t > 8)
   %  restart = t;
   %  break
   %end
   
   
 
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
	
