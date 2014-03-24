function g=trf(varargin)
% function g=trf(tw,[theta])
% Returns the (4 x 4) SE(3) rigid tranformation.
% The same as exp(tw,[theta])

% Kjartan Halvorsen
% 1999-06-01

tw=varargin{1};

if (nargin==1)
  theta=tw.theta;
else
  theta=varargin{2};
end

vw=tw.coordinates*theta;
v=vw(1:3);
v(:);
w=vw(4:6);

w_hat=vect2mat(w);
t=[w_hat v;0 0 0 0];

g=expm(t);