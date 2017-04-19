function [f, gradf] = gyro_bias_error(b, ws, q0, q1, dt)
  N = size(ws,2);
  qhats = zeros(N,4);
  phis = zeros(N,4);
  qhat = [0 0 0 1];
  for i=1:N
    phi_i = qexp( dt*(ws(:,i)-b) );
    qhat = quaternion_mult(qhat, phi_i);
    qhats(i,:) = qhat;
    phis(i,:) = phi_i';
  end
  q0qhat = qmult(q0, qhat);

  %%f = norm(qhat(:) - q1(:));
  diffq = qhat(:) - q1(:);
  f = 0.5*sum(diffq.^2);

  if (nargout == 2)
    %% Compute gradient
    qhatsR = zeros(N,4);
    qR = [0 0 0 1];
    for i=N:-1:1
      qR = qmult(phis(i,:), qR); 
      qhatsR(i,:) = qR;
    end

    dfdx = diffq;
    dphidb = [-dt/2*eye(3) zeros(3,1)]; % Approximation

    %% Special treatment of first and last element of the sum
    dqhatdb = zeros(3,4);
    for j=1:3
      dqhatdb(j,:) = dqhatdb(j,:) + qmult(dphidb(j,:), qhatsR(1,:));
    end
    for j=1:3
      dqhatdb(j,:) = dqhatdb(j,:) + qmult(qhats(N,:), dphidb(j,:));
    end

    for i=2:N-1
      for j=1:3
	dqhatdb(j,:) = dqhatdb(j,:) \
	    + qmult( qmult(qhats(i,:), dphidb(j,:)),  qhatsR(i,:));
      end
    end

    for j=1:3
      dqhatdb(j,:) = qmult(q0, dqhatdb(j,:));
    end
    
    gradf = (dqhatdb * dfdx);
  end
end
