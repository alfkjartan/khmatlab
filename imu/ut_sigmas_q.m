%UT_SIGMAS_q - Generate Sigma Points for Unscented Transformation for
%state vector with quaternion
%
% Syntax:
%   X = ut_sigmas(M,P,c);
%
% In:
%   M - Initial state mean (Nx1 column vector)
%   P - Initial state covariance
%   c - Parameter returned by UT_WEIGHTS
%
% Out:
%   X - Matrix where 2N+1 sigma points are as columns
%
% Description:
%   Generates sigma points and associated weights for Gaussian
%   initial distribution N(M,P). For default values of parameters
%   alpha, beta and kappa see UT_WEIGHTS.
%
% See also UT_WEIGHTS UT_TRANSFORM UT_SIGMAS
%

% Copyright (C) 2006 Simo Särkkä
%
% $Id: ut_sigmas.m 109 2007-09-04 08:32:58Z jmjharti $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

%% M contains [d, v, a, q, w, alpha, ed, ev, ea, eq, ew, eaa]

function X = ut_sigmas_q(M,P,c);

%  A = schol(P);
  A = chol(P)';
  n = size(M,1);
  X = zeros(n, 2*n+1);
  %X = sqrt(c)*[zeros(size(M)) A -A];
  AA = sqrt(c)*[A -A];
  MM = repmat(M,1,size(X,2));
  X(1:9,:) = MM(1:9,:) + cat(2, zeros(9,1), AA(1:9,:));
  X(10:13,2:end) = qpropagate(MM(10:13,2:end), AA(10:12,:), 1); 
  keyboard
  X(14:end, :) = MM(14:end,:) + cat(2, zeros(6,1), AA(13:end,:));

