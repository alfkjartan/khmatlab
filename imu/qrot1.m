function [vv, dvvdq] = qrot1(v, q, qc)

  if (nargin == 0)
    do_unit_test();
  else
    qv = qmult(q,v);
    vqc = qmult(v,qc);

    vv = qmult(qv, qc)';

    if (nargin > 1)
      I4 = eye(4);
      Ic4 = diag([-1 -1 -1 1]);
      
      dvvdq = zeros(4,4);

      for i=1:4
	dvvdq(:,i) = qmult(I4(:,i), vqc) + qmult(qv, Ic4(:,i));
      end
    end
  end

function do_unit_test()

  q = (quaternion(randn(3,1), rand(1)))';
  qc = qconj(q)';

  I4 = eye(4);
  dq = 1e-6;
  tol = 1e-5;

  v = cat(1, randn(3,1), 0);

  [vv, dvvdq] = qrot1(v, q, qc);
  v2 = qrot1(vv, qc, q);
  
  if ( norm(v-v2) > tol )
    disp('Test 1. Failed')
    cat(2, v, v2)
    keyboard
  else
    disp('Test 1. OK')
  end
  

  dvv_dq = zeros(4,4);

  for i=1:4
    dvv_dq(:,i) = (qrot1(v, q + I4(:,i)*dq, qconj(q + I4(:,i)*dq)) - vv) / dq;
  end

  if ( norm(dvvdq - dvv_dq) > tol )
    disp('Test 2. Failed')
    cat(2, dvvdq, dvv_dq)
    keyboard
  else
    disp('Test 2. OK')
  end

  
    
  
