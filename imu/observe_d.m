function [y]=observe_d(x)
%  [y]=observe_d(x, qi)
% Measurement function for an imu with orientation represented by
% a quaternion. 
%
% Input
%    x      ->   the current state. (n x N) vector.  x=[d,v,acc]

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

  y = x(7:9,:);
end

function do_unit_test()
