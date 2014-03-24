%UKF_UPDATEQ -  Unscented Kalman Filter update step for state vector
% with quaternion representation of rotation
%
% Syntax:
%   [M,P,K,MU,IS,LH] = UKF_UPDATE2(M,P,Y,h,R,qindstart,h_param,alpha,beta,kappa,mat)
%
% In:
%   M  - Mean state estimate after prediction step
%   P  - State covariance after prediction step
%   Y  - Measurement vector.
%   R  - Measurement covariance.
%%   qindstart      ->   Row index where quaternion starts. Empty matrix
%   param - Parameters of h               (optional, default empty)
%
% Out:
%   M  - Updated state mean
%   P  - Updated state covariance
%   K  - Computed Kalman gain
%   MU - Predictive mean of Y
%   S  - Predictive covariance Y
%   LH - Predictive probability (likelihood) of measurement.
%   
% Description:
%   Perform augmented form Discrete Unscented Kalman Filter (UKF)
%   measurement update step. Assumes additive measurement
%   noise.
%
%% Based on UKF_UPDATE2 in ekfukf-library
%
% See also:
%   UKF_PREDICT1, UKF_UPDATE1, UKF_PREDICT2, UKF_PREDICT3, UKF_UPDATE3
%   UT_TRANSFORM, UT_WEIGHTS, UT_MWEIGHTS, UT_SIGMAS

% History:
%   08.02.2008 JH Fixed a typo in the syntax description. 
%   04.05.2007 JH Initial version. Modified from ukf_update1.m
%              originally created by SS.
%   
% 
% References:
%   [1] Wan, Merwe: The Unscented Kalman Filter
%
% Copyright (C) 2007 Jouni Hartikainen, Simo S�rkk�
%
% $Id: ukf_update2.m 480 2010-10-18 07:45:48Z jmjharti $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M1,P1,v,K,Pxz, Pz, MU,S,LH] = ukf_updateq(M,P,Y,R,qi, p0, d0 , \
				       alpha,beta,kappa, mat)

  %
  % Check that all arguments are there
  %
  if nargin < 6
    error('Too few arguments');
  end
  if nargin < 7
    alpha = [];
  end
  if nargin < 8
    beta = [];
  end
  if nargin < 9
    kappa = [];
  end
  if nargin < 10
    mat = [];
  end

  %
  % Apply defaults
  %
  if isempty(mat)
    mat = 0;
  end

  %
  % Do transform and make the update
  %
  m = size(M,1); 
  n = size(R,1); 

  %% create sigma points
  [X,Wm,Wc] = sigmaq(P,M,qi);
  %% propagate sigma points
  %Z = observe_rb(X, p0, d0,qi);
  %Z = observe_rbb(X, p0, d0,qi);
  Z = observe_q(X, p0, d0,qi);


  [Xm, Zm, Px, Pz, Pxz] = sigma_cov_q(X, Z, qi, [], Wm, Wc);

  Pv = Pz + R;
  K = Pxz / Pv;
  v = Y-Zm;

  %% Update. Quaternion part is updated using quaternion multiplication
  M1 = qupdate(M, K, v, qi); 

  %xv = qstate2vstate(M, qi);
  %xvup = xv + K * v;
  %M1 = vstate2qstate(xvup, qi);
  

  P1 = P - K * Pv * K'; % P - Pxz* Pv^-1 * Pv * Pv^-1 * Pzx 
		        % = P - Pxz * Pv^-1 * Pzx

  try
    chol(P1);
  catch
    keyboard
  end

