function [d, vel, accs, q, ws, alphas] = track_imu(acc, w, dt, bandwidth)
% [d, vel, q] = track_imu(acc, w, dt, bandwidth)
% Tracks an imu from acceleration and gyro data.
% Uses part of the implementation of ukf in the EKF/UKF toolbox from Aalto university.
%
% Input
%    acc        ->   acceleration data (3 x nfr)
%    w          ->   gyro data (3 x nfr)
%    dt         ->   sampling period
%    bandwidth  ->   Filter tuning parameter
% Output
%    d          <-   the estimated translation of the rigid body.
%                    Represents translations from pos in first frame
%    vel        <-   velocity of the RB.
%    q          <-   quaternion representing the rotation of the RB.
%                    Rotates a point in the first frame to subsequent positions.

%% Kjartan Halvorsen
%% 2012-03-27. 
%%

%% Implemented as two sequential unscented Kalman filters
%%
%%   acc
%%  ----------------|                 ws   |-------------------------->
%%                  |----->|-------|-------|      
%%   w    |------|   q,w   | rot q |  accs       |------|  vs, ds
%%  ----->| QUKF |-------->|-------|------------>| DUKF |------------>
%%        |------|                               |------|
  if (nargin == 0)
    do_unit_test();
    %%do_unit_test();
  else

    n1 = 7;
    n2 = 9;

    %% The first state vector is x1=[q,w,alpha], or
    %%  x1=[q,w]
    %% The second state vector is x2=[d,v,acc]

    %%q1 = cat(2, [1 1 1]*pi/180*10, [1 1 1]*10*pi/180, [1 1 1]).^2;
    q1 = cat(2, [1 1 1]*pi/180*10, [1 1 1]*10*pi/180);
    qindstart = 1;
    q2 = cat(2, [0.01 0.01 0.01],[0.1 0.1 0.1],[1 1 1]);
    Q1 = diag(q1);
    Q2 = diag(q2);
  
    nfr=size(acc,2);

    sd_acc = 0.05 * 9.82; % in m/s^2
    sd_w = 10 * pi /180 ; % in rad/s
    R1 = diag([1 1 1]*sd_w.^2);
    R2 = diag([1 1 1]*sd_acc.^2);


    P1 = Q1*4;
    P2 = Q2*4;
    %%x1 = zeros(10, 1);
    x1 = zeros(7, 1);
    x1(5:7) = w(:,1);
    x1(qindstart:qindstart+3) = [0;0;0;1];

    x2 = zeros(9,1);
    x2(7:9) = acc(:,1);

    
    x2_pred = zeros(n1,nfr);
    x2_upd = zeros(n1,nfr);
    x2_pred = zeros(n2,nfr);
    x2_upd = zeros(n2,nfr);
    ws = zeros(3,nfr);

%  try
    for fr = 1:nfr
      [x1, P1] = ukf_predictq(x1, P1, Q1, qindstart, dt);
      x1_pred(:,fr) = x1;
      [x1, P1] = ukf_updateq(x1,P1,w(:,fr), R1, qindstart);
      x1_upd(:,fr) = x1;
      %qi = qinv(x1(1:4));
      qi = (x1(1:4));
      accs = qtransv(acc(:,fr), qi);
      ws(:,fr) = qtransv(w(:,fr), qi);
      [x2, P2] = ukf_predictd(x2, P2, Q2, dt);
      x2_pred(:,fr) = x2;
      [x2, P2] = ukf_update_d(x2,P2,accs, R2);
      x2_upd(:,fr) = x2;
    end
 % catch 
 %   keyboard
 % end
  

  d = x2_upd(1:3,:);
  vel = x2_upd(4:6,:);
  accs = x2_upd(7:9,:);
  q = x1_upd(1:4,:);
  %%alphas = x1_upd(8:10,:);

end

function do_unit_test()
  disp("Unit test for function track_imu")

  %% Create data of rigid body rotating about a fixed axis going through
  %% the origin and in the z-direction.
  %% Let the imu be located at [1;0;0] initially.

  N = 200;
  dt = 0.01;
  r = 0.5;
  q0 = [0;0;0;1];
  d0 = [r;0;0];
  v = [0;0;1];
  w0 = 200 * pi / 180; % rad/s
  %%ws = ( 1 + sin(linspace(0,4*pi, N)) )*w0;
  modf = (1 + sin(linspace(-pi/2,3*pi/2,N)));
  %wb = repmat([0;0;1]*w0, 1, N).*modf;
  wb = cat(1, zeros(2,N), w0*modf);
  ab = (centraldiff(wb',1/dt))';

  accs1 = zeros(3,N);

  for i=2:N
    wi = wb(:,i);
    accs1(:,i) = 0;
  end 
  %accb = repmat([-1;0;0]*r*w0.^2, 1, N);
  accb = cat(1, -wb(3,:).^2*r, cat(2, (wb(3,2:end)-wb(3,1:end-1))*r/dt, \
				   0), zeros(1,N));
  
  qs = zeros(4,N);
  w = zeros(3,N);
  d = zeros(3,N);
  
  qi = q0;
  for i=1:N
    qi = qpropagate(qi,wb(:,i),dt);
    qs(:,i) = qi;
    d(:,i) = -d0 + qtransv(d0,qi);
  end

  SNR = 10;
  eacc = randn(size(accb))'*r*w0^2/SNR;
  ew = randn(size(wb))'*w0/SNR;
  [b,a] = butter(4,0.8);
  eaccf = filtfilt(b,a,eacc);
  ewf = filtfilt(b,a,ew);

  [dd, vv, accs, qq, ws] = track_imu(accb+eaccf', wb+ewf', dt, 1);
  %[dd, vv, accs, qq, ws] = track_imu(accb+eaccf', wb, dt, 1);
  %[dd, vv, accs, qq, ws, alphas] = track_imu(accb, wb, dt, 1);

  dds = zeros(3,N);
  dds(:,1) = dd(:,1);

  errangle = zeros(N,1);
  accb2 = zeros(3,N);
  for i=1:N
    eq = qmult(qs(:,i), qinv(qq(:,i)));
    errangle(i) = 2*acos(eq(4))*180/pi;
    if (i>1)
      dds(:,i) = dds(:,i-1) + dt*vv(:,i-1);
    end 
    accb2(:,i) = qtransv(accs(:,i), qinv(qq(:,i)));
  end

  figure(1)
  clf
  plot(d(1,:), d(2,:));
  hold on
  plot(dd(1,:), dd(2,:), 'r');
  plot(dds(1,:), dds(2,:), 'm');

  d1 = [1;0;0];
  d1s = zeros(3,N);
  for i=1:N
    d1s(:,i) = qtransv(d1,qq(:,i));
  end

  figure(2)
  clf
  plot3(d1s(1,:), d1s(2,:), d1s(3,:))

  [ath, ar] = cart2pol(accs(1,:), accs(2,:));
  figure(4)
  clf
  subplot(211)
  plot(ath*180/pi)
  subplot(212)
  plot(ar)

  [vth, vr] = cart2pol(vv(1,:), vv(2,:));
  figure(5)
  clf
  subplot(211)
  plot(vth*180/pi)
  subplot(212)
  plot(vr)


  keyboard

function do_unit_test2()
  disp("Unit test for function track_imu")

  %% Create data of rigid body rotating about a fixed axis going through
  %% the origin and in the z-direction.
  %% Let the imu be located at [1;0;0] initially.

  N = 200;
  dt = 0.01;
  r = 0.5;
  q0 = [0;0;0;1];
  d0 = [r;0;0];
  v = [0;0;1];
  w0 = 200 * pi / 180; % rad/s
  %%ws = ( 1 + sin(linspace(0,4*pi, N)) )*w0;
  modf = (1 + sin(linspace(-pi/2,3*pi/2,N)));
  %wb = repmat([0;0;1]*w0, 1, N).*modf;
  wb = cat(1, zeros(2,N), w0*modf);
  

  %accb = repmat([-1;0;0]*r*w0.^2, 1, N);
  accb = cat(1, -wb(3,:).^2*r, zeros(2,N));
  
  qs = zeros(4,N);
  w = zeros(3,N);
  d = zeros(3,N);
  
  qi = q0;
  for i=1:N
    qi = qpropagate(qi,wb(:,i),dt);
    qs(:,i) = qi;
    d(:,i) = -d0 + qtransv(d0,qi);
  end

  SNR = 10;
  eacc = randn(size(accb))'*r*w0^2/SNR;
  ew = randn(size(wb))'*w0/SNR;
  [b,a] = butter(4,0.8);
  eaccf = filtfilt(b,a,eacc);
  ewf = filtfilt(b,a,ew);

  [dd, vv, accs, qq, ws] = track_imu(accb+eaccf', wb+ewf', dt, 1);
  %[dd, vv, accs, qq, ws] = track_imu(accb, wb+ewf', dt, 1);
  %[dd, vv, accs, qq, ws, alphas] = track_imu(accb, wb, dt, 1);

  dds = zeros(3,N);
  vvs = zeros(3,N);
  for i=2:N
    a1 = qtransv(accs(:,i-1), (qq(:,i-1)));
    a2 = qtransv(accs(:,i), (qq(:,i)));
    vvs(:,i) = vvs(:,i-1) + 0.5*dt* (a1 + a2); 
    dds(:,i) = dds(:,i-1) + 0.5*dt*(vvs(:,i-1) + vvs(:,i)) ...
	+ 0.25*dt^2*(a1 + a2);
  end

  figure(1)
  clf
  plot(d(1,:), d(2,:));
  hold on
  plot(dds(1,:), dds(2,:), 'r');

  d1 = [1;0;0];
  d1s = zeros(3,N);
  for i=1:N
    d1s(:,i) = qtransv(d1,qq(:,i));
  end

  figure(2)
  clf
  plot3(d1s(1,:), d1s(2,:), d1s(3,:))

  keyboard
