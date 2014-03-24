function xx = qupdate(x, K, v, nq)
  %% xx = qupdate(x, K, v, nq)
  %% Returns the updated (corrected) state estimate
  %%   xx = x + K*v
  %% Except for quaternion part which is
  %%   qq = q * qexp(Kq*vq)
  %%
  %% Input
  %%    x   -> predicted state estimate
  %%    K   -> Kalman filter gain
  %%    v   -> innovation (z - zm)
  %%    nq  -> Number of quaternions. Always in the first part of the
  %%           state vector
  %% Output
  %%    xx  <- updated state estimate

  %% Kjartan Halvorsen
  %% 2012-03-26

  tol = 1e-10;

  if (nq == 0)
    xx = x + K*v;
  else
    %% Contains quaternion
    xx = x;
    w = K*v; % The corrections to make
    xx(nq*4+1:end) = x(nq*4+1:end) + w(nq*3+1:end);
    for i=1:nq
      qix=(i-1)*4+1:i*4;
      qiw=(i-1)*3+1:i*3;
      xx(qix) = qmult(x(qix), qexp(w(qiw)*0.5));
    end
    
    %%Debug
    %if ( (norm(w(1:3)) < tol) & (norm(x(9:11)) > 1e-2))
    %  keyboard
    %end

end
  