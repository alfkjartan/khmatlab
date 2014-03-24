function ex=exp(varargin)
% function ex=exp(tw,theta)
% Computes the exponential of the twist.

% Kjartan Halvorsen
% 1999-05-31

thr=1e-12;

tw=varargin{1};

if (nargin==1)
  theta=tw.theta;
else
  theta=varargin{2};
end

t=gettwist(tw); % Returns the (4 x 4) twist
w_hat=t(1:3,1:3);
w=vect(w_hat);
v=t(1:3,4);

% The rotation
ew_hat=expm(t(1:3,1:3)*theta);

% The translation
if (norm(ew_hat-eye(3))<thr)
  d=v*theta;
else
  d=(eye(3)-ew_hat)*cross(w,v) + (w*w')*v*theta;
end

% Compiling
ex=[ew_hat d;zeros(1,3) 1];