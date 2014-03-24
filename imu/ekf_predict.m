%[x, S] =
%     ekf_predict(x0, S0, dyn_func, dyn_func_params, R)
%
% EKF using todorovs formulation and implementation. Made more general
% by using feval and name of  dynamics- and residual functions.
%
% x0: current state
% S0: covariance of current state estimate
%
% dyn_func, dyn_func_params: Name of function that computes the state
% prediction and its jacobian. 
%
% R: covariance of the dynamics noise: x(t+1) = x(t) + noise
%
%% Output
% x: predicted state
% S: covariance of the predicted state


function [x, S] = ekf_predict(x0, S0, dfunc, dfparams, R)

  [xx,A] = feval(dfunc, x0, dfparams);        % initialize estimate
   S = A*So*A' + R;                                % prior covariance
