function [d, vel, accs, q, ws, qf, df] = track_twonode(acc1, w1, acc2, w2, dt, bandwidth,q12,ekf)
% [d, vel, accs, q, ws, qf, df] = track_twonode(acc1, w1, acc2, w2, dt,
% bandwidth, ekf)
% Tracks a rigid body using data from two imus.
% Uses part of the implementation of ukf in the EKF/UKF toolbox from Aalto university.
%
% Input
%    acc1        ->   acceleration data (3 x nfr) from node 1
%    w1          ->   gyro data (3 x nfr) from node 1
%    acc2        ->   acceleration data (3 x nfr) from node 2
%    w2          ->   gyro data (3 x nfr) from node 2
%    dt          ->   sampling period
%    bandwidth   ->   Filter tuning parameter
%    ekf         ->   if not 0, use ekf. 
% Output
%    d          <-   the estimated translation of the rigid body.
%                    Represents translations from pos in first frame
%    accs       <-   acceleration in the spatial frame
%    vel        <-   velocity of the RB.
%    q          <-   quaternion representing the rotation of the RB.
%                    Rotates a point in the first frame to subsequent positions.
%    ws         <-   angular velocity in the spatial frame
%    qf         <-   orientation of node 2 wrt node 1
%    df         <-   position of node 2 wrt node 1

%% Kjartan Halvorsen
%% 2012-04-03. 
%%

%% Implemented as two sequential unscented Kalman filters
%%
%%   acc1
%%  --------------------|                 ws   |-------------------------->
%%                      |----->|-------|-------|      
%%   w1,w2    |------|   q,w   | rot q |  accs       |------|  vs, ds
%%  --------->| QUKF |-------->|-------|------------>| DUKF |------------>
%%            |------|                               |------|
  if (nargin == 0)
    do_unit_test();
    %%do_unit_test();
  else

    if (nargin < 8)
      ekf = 0;
    end

    if ekf
      resfunc1 = "residual_q2err";
      resfunc2 = "residual_d";
      Idt = eye(3)*dt;
      I3 = eye(3); Z3 = zeros(3,3);
      A1 = eye(6);
      A1(1:3,1:3) = 0;
      A2 = cat(1, \
	       cat(2, I3, Idt, Z3),\
	       cat(2, Z3, I3, Idt),\
	       cat(2, Z3, Z3, I3)); % Linear dynamics model for
				% translation part

      n1 = 6;
      qu1 = [0; 0; 0; 1];
      q12 = q12(:);
      qu12 = q12;
      x1 = zeros(n1, 1);
      x1(4:6) = w1(:,1);
      q1 = cat(2,  [5 5 5]*1e-1*pi/180, [1 1 1]*pi/180*10);

    else
      %% The first state vector is x1=[qf,q,w]
      n1 = 11;
      x1 = zeros(n1, 1);
      x1(1:4) = q12;
      x1(5:8) = [0;0;0;1];
      x1(9:11) = w1(:,1);
      q1 = cat(2,  [1 1 1]*1e-2*pi/180, [1 1 1]*pi/180*10, [1 1 1]*20*pi/180);
    end

    %% The second state vector is x2=[d,v,acc]
    n2 = 9;

    nq = 2;
    q2 = cat(2, [0.01 0.01 0.01],[0.1 0.1 0.1],[1 1 1]);

    Q1 = diag(q1.^2);
    Q2 = diag(q2.^2);
  
    nfr=size(acc1,2);

    sd_acc = 0.05 * 9.82; % in m/s^2
    %sd_w = 10 * pi /180 ; % in rad/s
    sd_w = 5 * pi /180 ; % in rad/s
    R1 = diag([1 1 1 1e1 1e1 1e1]*sd_w.^2);
    R2 = diag([1 1 1]*sd_acc.^2);


    P1 = Q1;
    P1(1:3,1:3) = eye(3)*(15/180*pi).^2;

    P2 = Q2*2;
    %%x1 = zeros(10, 1);

    x2 = zeros(9,1);
    x2(7:9) = acc1(:,1);

    
    x1_pred = zeros(n1,nfr);
    x1_upd = zeros(n1,nfr);
    x2_pred = zeros(n2,nfr);
    x2_upd = zeros(n2,nfr);
    ws = zeros(3,nfr);
    qf = zeros(4,nfr);
%  try
    pred1 = zeros(3,1);
    upd1 = zeros(3,1);
    pred2 = 0;
    upd2 = 0;
    for fr = 1:nfr
      if ekf
	x1(1:3) = 0; % The orientation of imu2 wrt imu1 considered constant. 
	P1 = A1*P1*A1' + Q1;
      else
	[x1, P1, slask, tt] = ukf_predictq2(x1, P1, Q1, nq, dt);
	pred1 = pred1 + tt;
      end
      %%keyboard
      x1_pred(:,fr) = x1;

      if ekf
	[x1,P1] = ekf_update(cat(1, w1(:,fr), w2(:,fr)), x1, P1, \
			     resfunc1, qu12, R1);
	qu1 = qpropagate(qu1, x1(4:6), dt);
	qu12 = qpropagate(qu12, x1(1:3), 1);
	q(:,fr) = qu1;
	qf(:,fr) = qu12;
	qi = qu1;
      else
	[x1, P1,tt] = ukf_updateq2(x1,P1, cat(1, w1(:,fr), w2(:,fr)), R1, nq);
	upd1 = upd1 + tt;
	qi = x1(5:8);
	qf(:,fr) = x1(1:4);
      end
			%keyboard
      x1_upd(:,fr) = x1;
      accs = qtransv(acc1(:,fr), qi);
      ws(:,fr) = qtransv(w1(:,fr), qi);

      tic();
      if ekf
	[x2] = dynamics_d(x2, dt);
	P2 = A2*P2*A2' + Q2;
      else
	[x2, P2] = ukf_predictd(x2, P2, Q2, dt);
      end

      pred2 = pred2 + toc();
      x2_pred(:,fr) = x2;

      tic();
      if ekf
	[x2,P2] = ekf_update(accs, x2, P2, \
			     resfunc2, [], R2);
      else
	[x2, P2] = ukf_update_d(x2,P2,accs, R2);
      end
      upd2 = upd2 + toc();
      x2_upd(:,fr) = x2;
    end
% catch 
%   keyboard
% end
  

  d = x2_upd(1:3,:);
  vel = x2_upd(4:6,:);
  accs = x2_upd(7:9,:);
  if ~ekf
    q = x1_upd(5:8,:);
  end

  %%alphas = x1_upd(8:10,:);

  %disp('Prediction 1: '), disp(pred1)
  %disp('Update 1: '), disp(upd1)
  %disp('Prediction 2: '), disp(pred2)
  %disp('Update 2: '), disp(upd2)
end

function   [wb1, wb2, accb1, accb2] = generate_data(N, w0, dt, axs, r, \
					   q12, d12, g);
  modf = 1 + sin(linspace(-pi/2,3*pi/2,N));
  wb = w0*modf;
  ab = (centraldiff(wb',1/dt))';

  wb1 = zeros(3, N);
  wb1(axs,:) = wb;
  wb2 = qrot(wb1,q12);
  
  zN = zeros(1,N);
  switch axs
    case 1
      accb1 =  cat(1, zN, -wb.^2*r, ab*r);
    case 2
      accb1 =  cat(1, ab*r, zN, -wb.^2*r);
    case 3
      accb1 =  cat(1, -wb.^2*r, ab*r, zN);
  end 
  accb2 = qrot(accb1,q12);

function do_unit_test()
  disp("Unit test for function track_twonode")
  debug = 0;
  ekf = 1;

  %% Create data of rigid body rotating about a fixed axis going through
  %% the origin and in the z-direction.
  %% Let the first imu be located at [1;0;0] initially and the second at
  %% [2;0;0]. The second is rotated about its x-axis 

  cycles = 6
  N = 100;
  w0 = 200 * pi / 180; % rad/s
  axes = [3,2];
  dt = 0.01;
  r = 0.5;
  g = 9.82*randn(3,1); %% The g vector
  d12 = 1;
  q12 = quaternion([1;0;0], pi/3);

  [wb1z, wb2z, accb1z, accb2z] = generate_data(N, w0, dt, 3, r, \
					   q12, d12, g);
  [wb1y, wb2y, accb1y, accb2y] = generate_data(N, w0, dt, 2, r, \
					   q12, d12, g);

  wb1 = repmat(cat(2, wb1z, wb1y), 1, cycles);
  wb2 = repmat(cat(2, wb2z, wb2y), 1, cycles);
  accb1 = repmat(cat(2, accb1z, accb1y), 1, cycles);
  accb2 = repmat(cat(2, accb2z, accb2y), 1, cycles);
  


  q12e = quaternion(randn(3,1),pi*10/180); % Accuracy in apriori angle measurement
  q0 = [0;0;0;1];
  d0 = [r;0;0];
  v = [0;0;1];


  N = size(wb1,2);

  qs = zeros(4,N);
  w = zeros(3,N);
  d = zeros(3,N);
  
  qi = q0;
  for i=1:N
    qi = qpropagate(qi,wb1(:,i),dt);
    qs(:,i) = qi;
    d(:,i) = -d0 + qtransv(d0,qi);
  end

  SNR = 10;
  eacc = randn(size(accb1))'*r*w0^2/SNR;
  eacc2 = randn(size(accb2))'*r*w0^2/SNR;
  ew = randn(size(wb1))'*5/180*pi;
  ew2 = randn(size(wb2))'*5/180*pi;
  [b,a] = butter(4,0.8);
  eaccf = filtfilt(b,a,eacc);
  eaccf2 = filtfilt(b,a,eacc2);
  ewf = filtfilt(b,a,ew);
  ewf2 = filtfilt(b,a,ew2);

  tic();
  [dd, vv, accs, qq, ws, qf, df] = track_twonode(accb1+eaccf', wb1+ewf', \
					 accb2+eaccf2', wb2+ewf2', dt, \
						 1, qmult(q12e,q12), ekf);
  %%[dd, vv, accs, qq, ws, qf, df] = track_twonode(accb1+eaccf', wb1, \	
  %%				 accb2+eaccf2', wb2, dt, \
  %%						 1, qmult(q12e,q12), ekf);


  toc()

  dds = zeros(3,N);
  dds(:,1) = dd(:,1);

  errangle = zeros(N,1);
  erranglef = zeros(N,1);
  accb2 = zeros(3,N);
  q21 = qinv(q12);
  for i=1:N
    eq = qmult(qs(:,i), qinv(qq(:,i)));
    eqf = qmult(qf(:,i), q21);
    errangle(i) = 2*acos(eq(4))*180/pi;
    erranglef(i) = 2*acos(eqf(4))*180/pi;
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
  %plot(dds(1,:), dds(2,:), 'm');
  legend('True path', 'Estimated path')

  d1 = [1;0;0];
  d1s = zeros(3,N);
  for i=1:N
    d1s(:,i) = qtransv(d1,qq(:,i));
  end

  figure(2)
  clf
  plot3(d1s(1,:), d1s(2,:), d1s(3,:))

  figure(3)
  clf
  tid = (1:N)*dt;
  plot(tid,errangle, 'b')
  hold on
  plot(tid, erranglef,'r')
  legend('orientation error', 'q12 error')
  ylabel('Degrees')
  xlabel('Time')

  if debug
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
  end

 keyboard

