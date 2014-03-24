function xx = qupdate(x, K, v, qi)
  %% xx = qupdate(x, K, v, qi)
  %% Returns the updated (corrected) state estimate
  %%   xx = x + K*v
  %% Except for quaternion part which is
  %%   qq = q * qexp(Kq*vq)
  %%
  %% Input
  %%    x   -> predicted state estimate
  %%    K   -> Kalman filter gain
  %%    v   -> innovation (z - zm)
  %%    qi  -> index into x where quaternion part begins
  %% Output
  %%    xx  <- updated state estimate

  %% Kjartan Halvorsen
  %% 2012-03-26

  if (isempty(qi))
    xx = x + K*v;
  else
    %% Contains quaternion
    xx = x;
    w = K*v; % The corrections to make
    xx(1:qi-1) = x(1:qi-1) + w(1:qi-1);
    xx(qi:qi+3) = qmult(x(qi:qi+3), qexp(w(qi:qi+2)*0.5));
    xx(qi+4:end) = x(qi+4:end) + w(qi+3:end);
end
  