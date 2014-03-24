function [r, pred, JJ] = residual_q2err(y, x, q12)
%% [r, pred, JJ] = residual_q2err(y, x)
%% Function that calculates the residual model for a two-gyro imu. The
%% measurements are angular velocities, so no problem with taking simple
%% vector differences: r = y - h(x).
% Input
%    y      ->   the current measurement (6 x 1)
%    x      ->   the current state. x = [v w], (6 x 1)

% Output
%    r      <-   the residual r = y - h(x)
%    pred   <-   the predicted model output h(x)
%    JJ     <-   the Jacobian dr/dx

% Kjartan Halvorsen
% 2012-05-03
%
%

if (nargin == 0)
  do_unit_test();
else
  
  [ev,dexpv_dv] = qexp(x(1:3));
  q12e = qmult(q12,ev);
  %%ev = [0;0;0;1]; % Since x(1:3) will always be zero after prediction step
  %%q12e = q12;

  
  pred = cat(1, x(4:6),  myqtransv(x(4:6), q12e));
  r = y - pred;

  if (nargout > 2)
				% Compute the Jacobian 
    dr1dv = zeros(3,3);
    dr1dw = -eye(3);
    v = x(1:3);

    wv = cat(1, x(4:6), 0);
    dr2dv = zeros(3,3);
    evc = qconj(ev);
    wvexpvc = qmult(wv,evc);
    expvwv = qmult(ev,wv);
    for i=1:3
      dr2dvi = - qtrans( qmult(dexpv_dv(:,i), wvexpvc) \
			+ qmult(expvwv, qconj(dexpv_dv(:,i))), q12);
      dr2dv(:,i) = dr2dvi(1:3);
    end

    dwv_dw = eye(3);
    dr2dw = zeros(3,3);
    for i=1:3
      dr2dw(:,i) = - qtransv( dwv_dw(:,i), q12e);
    end 

    JJ = [dr1dv dr1dw
	  dr2dv dr2dw];
  end % nargout > 2
end % if nargin == 0

function vv = myqtransv(v,q)
  vq = cat(1, v, 0);
  vvq = qmult(qmult(q,vq), qconj(q));
  vvq = vvq(:);
  vv = vvq(1:3);


function do_unit_test()
  disp("Unit test of function residual_q2err")
  tol = 1e-12;

  
  q12 = quaternion([0;0;1], pi/2)';
  v0 = [1,1,1]*4/180*pi;
  w0 = [0;0;1]*(-pi/2);
  v = zeros(3,1);
  y1 = w0;
  y2 = qtransv(w0, qmult(q12,qexp(0.5*v0)));
  y = cat(1, y1, y2);
  x = cat(1, v, w0);
  [r, yy, JJ] = residual_q2err(y, x, q12);
  
  dx = eye(6);
  dw = 1e-6;
  drdw = zeros(6,3);
  for i=4:6
    drdw(:,i-3) = (residual_q2err(y, x+dx(:,i)*dw, q12) - r) / dw;
  end

  if ( norm(JJ(:,4:6) - drdw) > dw*2 )
    disp('Test 1. Failed')
    cat(2, JJ(:,4:6), drdw)
    keyboard
  else
    disp('Test 1. OK')
  end

  drdv = zeros(6,3);
  for i=1:3
    drdv(:,i) = (residual_q2err(y, x+dx(:,i)*dw, q12) - r) / dw;
  end

  if ( norm(JJ(:,1:3) - drdv) > dw*2 )
    disp('Test 2. Failed')
    cat(2, JJ(:,1:3), drdv)
    keyboard
  else
    disp('Test 2. OK')
  end
