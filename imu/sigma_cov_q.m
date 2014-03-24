function [mx, my, Sx, Sy, Cxy, EE] = sigma_cov_q(XX, YY, qix, qiy, Wm, Wc)
%%  [mx, my, Sx ,Sy, Cxy, EE] = sigma_cov_q(XX, YY, qindstartx, qindstarty, Wm, Wc)
%% Computes mean and covariance of sigma points for the system with quaternion representation
%% of rotation in elements 10:13 of the state vector.
%%
%% Input
%%   XX           ->   Sigma points 
%%   YY           ->   Sigma points. Could be empty 
%%   qindstartx    ->   Row index into XX where quaternion starts 
%%   qindstarty    ->   Row index into YY where quaternion starts 
%%
%% Output
%%   mu     <-   The mean of the sigma points 
%%   Sx     <-   The covariance of the set XX 
%%   Sy     <-   The covariance of YY
%%   Cxy      <-   The cross-covariance 



%% Compute the mean and separate covariances

if (nargin ==0)
  do_unit_test();
else

  [mx, Sx, EEx] = sigma_mean_q(XX,qix,Wm,Wc);
  [my, Sy, EEy] = sigma_mean_q(YY,qiy,Wm,Wc);

  nx = size(Sx,1);
  ny = size(Sy,1);

  dx = size(mx,1);
  dy = size(my,1);

  Cxy = zeros(nx,ny);

  %% Compute cross covariance

  for i=1:size(XX,2)
    if (isempty(qix) == 0) % State vector has quaternion
      Wx = cat( 1, XX(1:qix-1,i) - mx(1:qix-1), EEx(:,i), XX(qix+4:dx,i) - mx(qix+4:dx) );
    else
      Wx = XX(:,i) - mx;
    end
    if (isempty(qiy) == 0) % State vector has quaternion
      Wy = cat( 1, YY(1:(qiy-1),i) - my(1:qiy-1), EEy(:,i), YY(qiy+4:dy,i) - my(qiy+4:dy) );
    else
      Wy = YY(:,i) - my;
    end

    Cxy = Cxy + Wc(i) * Wx * Wy';
  end
end

function do_unit_test()
  tol = 1e-12;

  disp("Unit test of function sigma_cov_q")
  
  %% Generate some sigma points
  qi = 7;
  X1 = randn(13,1);
  X1(qi:qi+3) = quaternion(randn(3,1), rand(1));
  Sx = zeros(12,12);
  Sx(1,1) = 1;
  Sx(4,4) = 1;

  Xd = cat(1, Sx, zeros(1,12));
  X = repmat(X1, 1, 25) + cat(2, zeros(13,1), Xd, -Xd);

  p0 = randn(3,4);
  Y = observe_rb(X, p0, qi);

  [Wm,Wc] = ut_weights(12);

  [mx,my,Sx,Sy,Cxy] = sigma_cov_q(X,Y,qi,[], Wm, Wc);

  CCxy = Cxy;
  CCxy(1,[1:3:end]) = 0;
  if (norm(CCxy) > tol)
    disp("Test 1: OK")
    disp(sprintf("Unexprected result. Norm = %0.8f", norm(Cxy)))
    Cxy
  else
    disp("Test 1: OK")
  end

  X0 = zeros(13,1);
  X0(qi+3) = 1;
  X = sigmaq(eye(12), X0, qi);
  X(qi:qi+2,:) = 0;
  X(qi+3,:) = 1;
  Y = observe_rb(X, p0, qi);
  [mx,my,Sx,Sy,Cxy] = sigma_cov_q(X,Y,qi,[], Wm, Wc);
  p0myY=cat(2, p0(:), my, Y(:,1))
  
  grad = pi/180;
  X0 = zeros(13,1);
  X0(qi+3) = 1;
  X = sigmaq(eye(12)*5*grad, X0, qi);
  X([1:qi-1 qi+4:13],:) = 0;
  Y = observe_rb(X, p0, qi);
  [mx,my,Sx,Sy,Cxy] = sigma_cov_q(X,Y,qi,[], Wm, Wc);
  p0myY=cat(2, p0(:), my, Y(:,1))
  

  keyboard


    


