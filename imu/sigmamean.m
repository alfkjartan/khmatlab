function [mu,S, EEs, WW] = sigmamean(XX, nq, Wm, Wc)
%%  [mu,Sx, EEs] = sigmamean(XX, nq, Wm, Wc)
%% Computes mean and covariance of sigma points for the system with quaternion representation
%% of rotation in elements qindstart:qindstart+3 of the state vector.
%%
%% Input
%%   XX           ->   Sigma points 
%%   nq           ->   Number of quaternions in the state vector
%%
%% Output
%%   mu     <-   The mean of the sigma points 
%%   Sx     <-   The covariance of the set XX 
%%   EE     <-   Error vectors form calculation of covariance of
%%               quaternion part

%% Revisions

%% Compute the mean

  [n,nsigm] = size(XX);

  mu = zeros(n,1);
  for i=1:size(XX,2)
    mu = mu + Wm(i) * XX(:,i);
  end

  
  %% The quaternion part must be handled separately
  if (nq)  % State vector has quaternion
    EEs = zeros(3*nq,nsigm);
    for j=1:nq
      qix = (j-1)*4+1:j*4;
      qiw = (j-1)*3+1:j*3;
      [mu(qix), EEs(qiw,:), its] = qmean(XX(qix,:),Wm);
      %%disp('its'), disp(its)
    end
  else
    EEs = [];
  end

  %% and the covariance
  if (nq) % State vector has quaternion
    WW = zeros(n-nq, nsigm);
    WW(1:nq*3,:) = EEs;
    WW(nq*3+1:n-nq,:) = XX(nq*4+1:n,:) - repmat(mu(nq*4+1:n),1,nsigm);
    S = zeros(n-nq,n-nq);
    for i=1:nsigm
      W = WW(:,i);
      S = S + Wc(i) * W * W';
    end
  else
    S = zeros(n,n);
    WW = XX - repmat(mu,1,nsigm);
    for i=1:size(XX,2)
      W = XX(:,i) - mu;
      S = S + Wc(i) * W * W';
    end
  end


