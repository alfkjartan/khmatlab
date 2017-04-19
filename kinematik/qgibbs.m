function [q, dq_da, dqconj_da] = qgibbs(a)
%%  q = qgibbs(a)
%% Returns the quaternion corresponding to the gibbs vector a
%%

%% Kjartan Halvorsen
%% 2012-06-19

if (nargin == 0)
  do_unit_test();
else

  a2 = sum(a.^2);
  a2p4 = a2+4;
  nrm = 1/sqrt(4+a2);
  q = [a;2]*nrm;

  if (nargout > 1)
    f1 = nrm/a2p4*a';
    dq_da = [eye(3)*nrm;zeros(1,3)] ...
	- [a;2]*f1;
  end
  
  if (nargout > 2)
    dqconj_da = -dq_da - [zeros(3,3);4*f1];
  end

  
end

function do_unit_test()

  w = randn(3,1);
  dw = 1e-6;
  I3 = eye(3);

  tic()
  for i=1:1000
    [ew,JJ,JC] = qexp(w);
    [ewc,JJC] = qexp(-w);
  end
  toc()

  tic()
  for i=1:1000
    [ew,JJ,JC] = qgibbs(w);
    [ewc,JJC] = qgibbs(-w);
  end
  toc()
  
  de_dw = zeros(4,3);
  for i = 1:3
    de_dw(:,i) = (qgibbs(w+dw*I3(:,i)) - ew) / dw;
  end

  if ( norm(JJ - de_dw) > dw )
    disp('Test 1. Failed')
    cat(2, JJ, de_dw)
    keyboard
  else
    disp('Test 1. OK')
  end

  dec_dw = zeros(4,3);
  for i = 1:3
    dec_dw(:,i) = (qgibbs(-(w+dw*I3(:,i))) - ewc) / dw;
  end

  if ( norm(JC - dec_dw) > dw )
    disp('Test 2. Failed')
    cat(2, JC, dec_dw)
    keyboard
  else
    disp('Test 2. OK')
  end 

  [q, dqda, dqconjda] = qgibbs([0;0;0]);
  bv = [1;2;3;0];
  dgbda = zeros(4,3);
  for i=1:3
      dgbda(:,i) = (qmult(dqda(:,i), bv) + qmult(bv,dqconjda(:,i)))';
  end
  dgbda
  hat(bv(1:3))
  keyboard
end

endfunction

  
