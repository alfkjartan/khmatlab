%[X, Err, S, Res, Pred, Cost] = 
%     estimate(Y, x0, S0, segment, map, info, R, V)
%
% Motion estimation from sensor data, Gauss-Newton method
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
   estimate(Y, x0, S0, segment, map, info, R, V)

maxIter = 10;                                % max Gauss-Newton iter
sclFactor = 0.1;                             % scale factor for linesearch
sclMin = 1E-10;                              % stop linesearch if scl<min
minImprove = 0.01;                           % min relative improvement

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

   err = r(flagR)'*Vinv*r(flagR)/2;          % initial cost

   %----------- Gauss-Newton method with backtracking linesearch ----
   for it = 1:maxIter,                       % until maxIter or converge
      JV = J(flagR,:)'*Vinv;
      grad = Pinv*(x-x0);                    % compute gradient
      grad(flagJ)= grad(flagJ)+ JV*r(flagR);
      H = Pinv;                              % compute Hessian
      H(flagJ,flagJ) = H(flagJ,flagJ) + JV*J(flagR,:);
      Hinv = inv(H);                         % invert Hessian

      update = - Hinv*grad;                  % compute full update

      scl = 1;
      while scl>=sclMin,                     % backtracking linesearch
         newx = x + scl*update;              % update at current scale
         newr = residual(Y(:,d), newx, segment, map);
         newerr = newr(flagR)'*Vinv*newr(flagR)/2 + ...
            (newx-x0)'*Pinv*(newx-x0)/2;

         if newerr < err,
            break;                           % exit if error decreased
         else
            scl = scl*sclFactor;             % otherwise reduce scale
         end
      end

      derr = err - newerr;
      if derr > 0,                           % improvement found
         x = newx;                           % update solution
         err = newerr;
         [r, predict, J] = residual(Y(:,d), x, segment, map);
      end
      
      if derr < err*minImprove,              % improvement too small
         break;                              % exit
      end
   end

   %------------- complete the update ---------------------------------
   H = Pinv;                                 % recompute Hessian
   H(flagJ,flagJ) = H(flagJ,flagJ) + J(flagR,:)'*Vinv*J(flagR,:);
   S = inv(H);                               % posterior covariance

   X(:,d) = x;                               % assign estimate
   Err(:,d) = sqrt(diag(S));                 % assign standard error
   Res(:,d) = r;                             % assign residual
   Pred(:,d) = predict;                      % assign prediction
   Cost(d) = err;                            % assign cost

   x0 = x;                                   % prepare for next iteration
end
fprintf('\n');
