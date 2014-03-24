function [r, pred, JJ] = residual_d(y, x)
%% [r, pred, JJ] = residual_d(y, x)
%% Function that calculates the residual model for the position of an imu.
%% measurements are acceleration, so no problem with taking simple
%% vector differences: r = y - h(x).
% Input
%    y      ->   the current measurement (3 x 1)
%    x      ->   the current state. (9 x 1)

% Output
%    r      <-   the residual r = y - h(x)
%    pred   <-   the predicted model output h(x)
%    JJ     <-   the Jacobian dr/dx

% Kjartan Halvorsen
% 2012-04-24
%
%

[pred] = observe_d(x);
r = y - pred;
JJ = cat(2, zeros(3,6), -eye(3));

