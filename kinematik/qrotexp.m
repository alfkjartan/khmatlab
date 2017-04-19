function [vv, dvvdw] = qrotexp(v, q, w, qc, qe, qec, dqdw)

  if (nargin == 0)
    do_unit_test();
  else

    if (nargin <7)
      [ew, dewdw] = qexp(w);
      qe = qmult(q, ew);
      qec = qconj(qe);

      dqdw = zeros(4,3);
      for i=1:3
	dqdw(:,i) = qmult(q, dewdw(:,i));
      end

    end
    
    [vv, dvvdq] = qrot1(v, qe, qec);

    dvvdw = dvvdq*dqdw;
  end

function do_unit_test()

  q = (quaternion(randn(3,1), rand(1)))';
  qc = qconj(q)';
  w = randn(3,1);
  I3 = eye(3);
  dw = 1e-6;
  tol = 1e-5;

  v = cat(1, randn(3,1), 0);

  [vv, dvvdw] = qrotexp(v, q, w, qc);

  dvv_dw = zeros(4,3);

  for i=1:3
    dvv_dw(:,i) = (qrotexp(v, q, w+I3(:,i)*dw, qc) - vv) / dw;
  end

  if ( norm(dvvdw - dvv_dw) > tol )
    disp('Test 1. Failed')
    cat(2, dvvdw, dvv_dw)
    keyboard
  else
    disp('Test 1. OK')
  end

  
    
  
