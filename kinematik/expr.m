function [ex, d] = expr(M, theta)
% function [ex, d] = expr(M, theta)
% Computes the matrix exponential of the matrix M using Rodrigues'
% formula. 
% Based on twist/exp.m
% If two output arguments are required, then ex is the rotation part
% (3x3) and d the translation.
% If only one output argument, this is the 4x4 transformation matrix.
  
% Kjartan Halvorsen
% 2002-06-24

thr=1e-12;

switch ( size(M, 1) )
  case 6
   w = M(4:6);
   v = M(1:3);
   w_hat = hat(w);
  case 4
    w_hat=M(1:3,1:3);
    vw = vee(M);
    w=vw(4:6);
    v=vw(1:3);
  case 3
    w_hat=M(1:3,1:3);
    w=vect(w_hat);
    v = zeros(3,1);
end

if nargin == 1 %% Otherwise it is assumed that w is unit norm!!
  theta=norm(w);
  if theta>thr
    w=w/theta;
    w_hat=w_hat/theta;
    v=v/theta;
  end
end


% The rotation
ew_hat=expm_rodrigues(w_hat,theta);

% The translation
if (norm(ew_hat-eye(3))<thr)
  d=zeros(3,1);
else
  d=(eye(3)-ew_hat)*cross(w,v) + (w*w')*v*theta;
end

% Finally

if (nargout == 1)
  if size(M,1) == 3
    ex = ew_hat;
  else
    ex=[ew_hat d;zeros(1,3) 1];
  end
else
  ex = ew_hat;
end
