%[X, Err, S, Res, Pred, Cost] =
%     estimate_ekf(Y, x0, S0, dyn_func, dyn_func_params. resid_func,
%     resid_func_params, R, V)
%
% EKF using todorovs formulation and implementation. Made more general
% by taking name of  dynamics- and residual functions.
%
% Y: matrix of measurements, each column is a data vector
% x0: mean of the initial state
% S0: covariance of the initial state
%
% dyn_func, dyn_func_params: Name of function that computes the state
% prediction and its jacobian. 
% resid_func, resid_func_params: Name of residual function and its jacobian. 
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
   estimate_ekf(Y, x0, S0, dfunc, dfparams, rfunc, rfparams, R, V)

nD = size(Y,2);                              % length of data sequence
nX = length(x0);
nR = size(V,1);

X = zeros(nX, nD);                       % allocate state estimates
Err = zeros(nX, nD);                     % allocate standard errors
Res = zeros(nR, nD);                     % allocate residuals
Pred = zeros(nR, nD);                    % allocate predictions
Cost = zeros(1, nD);                         % allocate costs


%--------------- loop over data sequence, estimate ----------------------
S = S0;
x = x0;
for d = 1:nD
   if mod(d,100)==0,                         % print every 100 samples
      fprintf('.');
   end
   

   [xx,A] = feval(dfunc, x0, dfparams);        % initialize estimate
   
   P = A*S*A' + R;                                % prior covariance
   Pinv = inv(P);                            % precompute P inverse

   [r, predict, J] = ...
      feval(rfunc, xx, rfparams);

   flagR = (~isnan(r));                      % find available residuals
   Vinv = inv(V(flagR,flagR));               % invert relevant part of V

   %----------- extended Kalman filter ---------------------------------
   JV = J(flagR,:)'*Vinv;
   %%grad = Pinv*(x-x0);                       % compute gradient
   %%grad(flagJ)= grad(flagJ)+ JV*r(flagR);
   grad = zeros(size(x));                       % compute gradient
   grad(flagJ)= grad(flagJ)+ JV*r(flagR);
   H = Pinv;                                 % compute Hessian
   H(flagJ,flagJ) = H(flagJ,flagJ) + JV*J(flagR,:);
   Hinv = inv(H);                            % invert Hessian

   x = xx - Hinv*grad;                        % compute update
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
