function [r, pred, JJ] = residual_q2(y, x)
%% [r, pred, JJ] = residual_q2(y, x)
%% Function that calculates the residual model for a two-gyro imu. The
%% measurements are angular velocities, so no problem with taking simple
%% vector differences: r = y - h(x).
% Input
%    y      ->   the current measurement (6 x 1)
%    x      ->   the current state. (11 x 1)

% Output
%    r      <-   the residual r = y - h(x)
%    pred   <-   the predicted model output h(x)
%    JJ     <-   the Jacobian dr/dx

% Kjartan Halvorsen
% 2012-04-24
%
%

[pred, dhdx] = observe_q2(x);
r = y - pred;
JJ = -dhdx;
