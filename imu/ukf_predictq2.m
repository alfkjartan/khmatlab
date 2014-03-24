% Based on 
%UKF_PREDICTq  Augmented (state and process noise) UKF prediction step
%
% Syntax:
%   [M,P] = UKF_PREDICTQ(M,P,Q,dt,qindstart)
%
% In:
%   M - 19x1 mean state estimate of previous step
%   P - 18x18 state covariance of previous step
%   Q - Non-singular 18x18 covariance of process noise w
%   dt - sampling time
%%   qindstart      ->   Row index where quaternion starts. Empty matrix
%%   means no quaternion
%
% Out:
%   M - Updated state mean
%   P - Updated state covariance
%
% Description:
%   Perform Unscented Kalman Filter prediction step
%   for model with quaternion representation of rotation
%

% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [XXm,Pk, EE, timeused] = ukf_predictq2(M,P,Q,nq,dt)


  timeused = zeros(3,1);
  %% Create sigma points
  n= size(P,1);
  [Wm,Wc,c] = ut_weights(n);

  tic();
  X = sigma(P+Q,M,nq);
  timeused(1) = toc(); % about 3.0 / 200 s

  %% Propagate sigma points
  tic();
  XX = dynamics_q2(X,dt);
  timeused(2) = toc(); % about 1.6 / 200 s

  tic();
  [XXm,Pk,EE] = sigmamean(XX,nq,Wm,Wc);
  timeused(3) = toc(); % About 3.2 / 200 s

  
  

  

