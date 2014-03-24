function [y]=observe_q(x)
%  [y]=observe_q(x, qi)
% Measurement function for an imu with orientation represented by
% a quaternion. 
%
% Input
%    x      ->   the current state. (n x N) vector.  x=[q,w,alpha]

% Output
%    y     <-   the measurements [acc;w]

% Kjartan Halvorsen
% 2012-03-27
%
%

if (nargin == 0)
  do_unit_test();
else
  N = size(x,2);

  %% The angular velocity are given in the body frame.
  
  y = x(5:7,:);

end

function do_unit_test()
