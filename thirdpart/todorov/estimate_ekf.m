%[X, Err, S, Res, Pred, Cost] =
%     estimate_ekf(Y, x0, S0, segment, map, info, R, V)
%
% Motion estimation from sensor data, EKF method
%
% Y: matrix of measurements, each column is a data vector
% x0: mean of the initial state
% S0: covariance of the initial state
%
% segment, map, info: descripion of the kinematic model;
%                     see residual.m and prepare.m for details
%
% R: covariance of the dynamics noise: x(t+1) = x(t) + noise
% V: covariance of the residuals (i.e. sensor noise)
%
% X: matrix of estimated states, each column is a state vector
% Err: matrix of standard errors, same size as X
% S: covariance of the final state
% Res: matrix of residuals
% Pred: matrix of predicted measurements, same size as Y
% Cost: vector of log-likelihood costs at state estimates

% Copyright (C) Emanuel Todorov, 2006-2007

function [X, Err, S, Res, Pred, Cost] = ...
   estimate_ekf(Y, x0, S0, segment, map, info, R, V)

nD = size(Y,2);                              % length of data sequence

X = zeros(map.nX, nD);                       % allocate state estimates
Err = zeros(map.nX, nD);                     % allocate standard errors
Res = zeros(map.nR, nD);                     % allocate residuals
Pred = zeros(map.nY, nD);                    % allocate predictions
Cost = zeros(1, nD);                         % allocate costs

if map.nJ == map.nP+map.nW,                  % full estimation
   flagJ = (info.type==1 | info.type==2);    % flag p,w-elements of x
else                                         % time-varying only
   flagJ = (info.type==1);                   % flag p-elements of x   
end

%--------------- loop over data sequence, estimate ----------------------
S = S0;
for d = 1:nD
   if mod(d,100)==0,                         % print every 100 samples
      fprintf('.');
   end
   
   P = S + R;                                % prior covariance
   Pinv = inv(P);                            % precompute P inverse

   x = x0;                                   % initialize estimate
   [r, predict, J] = ...
      residual(Y(:,d), x, segment, map);     % compute r, J at init.

   flagR = (~isnan(r));                      % find available residuals
   Vinv = inv(V(flagR,flagR));               % invert relevant part of V

   %----------- extended Kalman filter ---------------------------------
   JV = J(flagR,:)'*Vinv;
   grad = Pinv*(x-x0);                       % compute gradient
   grad(flagJ)= grad(flagJ)+ JV*r(flagR);
   H = Pinv;                                 % compute Hessian
   H(flagJ,flagJ) = H(flagJ,flagJ) + JV*J(flagR,:);
   Hinv = inv(H);                            % invert Hessian

   x = x - Hinv*grad;                        % compute update
   S = Hinv;                                 % posterior covariance

   X(:,d) = x;                               % assign estimate
   Err(:,d) = sqrt(diag(S));                 % assign standard error

   % evaluate at estimated state if necessary
   if nargout>3,
      [r, predict] = residual(Y(:,d), x, segment, map);
      Res(:,d) = r;
      Pred(:,d) = predict;
      Cost(d)= r(flagR)'*Vinv*r(flagR)/2 + (x-x0)'*Pinv*(x-x0)/2;
   end
   
   x0 = x;                                   % prepare for next iter.   
end
fprintf('\n');
