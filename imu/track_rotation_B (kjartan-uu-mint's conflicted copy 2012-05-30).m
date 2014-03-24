function [x,T, q, q12e, p12e] = track_rotation_B(time1, w1, acc1, time2, w2, \
					acc2, bandwidth, q12, p12, \
						 fixp12, fixq12)
%  [x,T, q12e] = track_rotation_B(time1, w1, acc1, time2, w2, \
%%					acc2, bandwidth, q12, p12)
% Tracks the rotation of a rigid body using data from two gyros and two accelerometers using
% algorithm B 
%
% Input
%    time1       ->   time stamps from node 1
%    w1          ->   gyro data (3 x nfr) from node 1
%    acc1        ->   accelerometer data (3 x nfr) from node 1
%    time1       ->   time stamps from node 1
%    w2          ->   gyro data (3 x nfr) from node 2
%    acc2        ->   accelerometer data (3 x nfr) from node 2
%    bandwidth   ->   Filter tuning parameter
%    q12         ->   Quaternion. Initial guess of orientation between
%                     the nodes
% Output
%    x           <-   the estimated state vectors (6 x nfrs)
%    T           <-   the time stamp of each estimate

%% Kjartan Halvorsen
%% 2012-05-21. 
%%

 if (nargin == 0)
   do_unit_test();
   %%do_unit_test();
 else

   if (nargin < 10)
     fixp12 = 0;
   end

   if (nargin < 11)
     fixq12 = 0;
   end

   
   dt0 = (time1(100) - time1(99))/1000; 
   dt0 = 0.01; % hardcoded.

   resfunc = "residual_twoimu";
   if fixp12
     Q = zeros(12,12);
     P = zeros(12.12);
     Q(10:12, 10:12) = dt0*eye(3)*(1*pi/180)^2;
     P(10:12,10:12) = eye(3)*(10*pi/180)^2; %% Covariance of initial guess of
     %% orientation between sensors
     resparams = {q12, p12};
     x0 = zeros(12,1);
     x0(1:3) = w1(:,1);
     x0(7:9) = acc1(:,1);
     A = eye(12,12);
     A(10:12,10:12) = 0;
   else
     Q = zeros(15,15);
     P = zeros(15,15);
     Q(10:12, 10:12) = dt0*eye(3)*(1*pi/180)^2;
     Q(13:15, 13:15) = dt0*eye(3)*(0.001)^2;
     P(10:12, 10:12) = eye(3)*(10*pi/180)^2; %% Covariance of initial guess of
     P(13:15,13:15) = eye(3)*(0.002)^2; %% Covariance of initial guess of
     resparams = {q12};
     x0 = zeros(15,1);
     x0(1:3) = w1(:,1);
     x0(7:9) = acc1(:,1);
     x0(13:15) = p12;
     A = eye(15,15);
     A(10:12,10:12) = 0;
   end
   if fixq12
     P = zeros(9,9);
     Q = zeros(9,9);
     resparams = {q12, p12, fixq12};
     x0 = zeros(9,1);
     x0(1:3) = w1(:,1);
     x0(7:9) = acc1(:,1);
     A = eye(9,9);
     A(1:3,4:6) = eye(3)*dt0;
   end

   %% Use cont white noise acc model
   I3a = eye(3)*(1000)^2; %% (grad/s/s)^2
   Q(1:3,1:3) = 0.3333*dt0^3*I3a;
   Q(1:3,4:6) = 0.5*dt0^2*I3a;
   Q(4:6,1:3) = 0.5*dt0^2*I3a;
   Q(4:6,4:6) = dt0*I3a;
   
   Q(7:9, 7:9) = dt0*eye(3)*(100^2); % Acc

   P(1:9, 1:9) = Q(1:9, 1:9)/dt0;
   P(4:6,4:6) = eye(3);

   nfr1 = length(time1);
   nfr2 = length(time2);

   sd_w = 20 * pi /180 ; % in rad/s
   sd_acc = 1 ; % in m/s/s
   R = zeros(12,12);
   R(1:3,1:3) = eye(3)*sd_w^2;
   R(4:6,4:6) = eye(3)*sd_acc^2;
   R(7:9, 7:9) = R(1:3,1:3);
   R(10:12, 10:12) = R(4:6,4:6);

   x = zeros(length(x0), nfr1+nfr2);
   T = nan(nfr1+nfr2,1);
   q12e = zeros(nfr1+nfr2,4);
   p12e = zeros(3,nfr1+nfr2);
   q = zeros(nfr1+nfr2, 4);
   qq = [0 0 0 1]';

   currentTime = min(time1(1), time2(1))-dt0;
   xc = x0; % current estimate
   fr1 = 1;
   fr2 = 1;
   nan6 = nan(6,1);


   k=0;
   while ( (fr1 <= length(time1)) & (fr2 <= length(time2)))
     t1 = time1(fr1);
     t2 = time2(fr2);
     if (t1 < t2)
       dt = t1 - currentTime;
       currentTime = t1;
       y = cat(1, w1(:,fr1), acc1(:,fr1), nan6);
       %%y = cat(1, nan(3,1), acc1(:,fr1), nan6); % Test with no gyro data
       fr1 = fr1 + 1;
     else 
       %% Same timestamp or next data from sensor 2 arrives first
       if (t1 == t2)
	 dt = t2 - currentTime;
	 currentTime = t2;
	 y = cat(1, w1(:,fr1), acc1(:, fr1), w2(:,fr2), acc2(:,fr2));
	 %%y = cat(1, nan(9,1), acc2(:,fr2)); % Test with no gyro data
	 fr2 = fr2 + 1;
	 fr1 = fr1 + 1;
       else
	 dt = t2 - currentTime;
	 currentTime = t2;
	 y = cat(1, nan6, w2(:,fr2), acc2(:,fr2));
	 fr2 = fr2 + 1;
	 %%y = cat(1, nan(9,1), acc2(:,fr2)); % Test with no gyro data
       end
     end

     for its = 1:ceil(dt/dt0)
       xc = A*xc;
       P = A*P*A' + Q;
     end
     %%keyboard
     resparams{1} = q12;
     [xc,P] = ekf_update(y, xc, P, \
			 resfunc, resparams, R);
     if (length(xc)>9)
       q12 = qmult(q12, qexp(xc(10:12)));
     end
     qq = qmult(qq, qexp(dt*xc(1:3)));
     k = k+1;
     x(:,k) = xc;
     T(k) =currentTime;
     q(k,:) = qq;
     q12e(k,:) = q12;
   end
   q = q';
   q12e = q12e';
   keep = find(~isnan(T));
   x = x(:, keep);
   T = T(keep);
   q = q(:,keep);
   q12e = q12e(:,keep);
 end 

function   [wb1, accb1] = generate_data(N, w0, dt, axs, r, \
					   q12, g);
  modf = 1 + sin(linspace(-pi/2,3*pi/2,N));
  wb = w0*modf;
  time1 = (0:N-1)*dt;
  angleb = cumtrapz(time1, wb);
  %% Scale wb so that it becomes a complete revolution
  wb = wb*2*pi/angleb(end);

  ab = (centraldiff(wb',1/dt))';

  wb1 = zeros(3, N);
  wb1(axs,:) = wb;
  
  ab1 = zeros(3,N);
  ab1(axs,:) = ab;

  if (length(r) == 3)
    
    accb1 = zeros(3, N);
    for i=1:N
      accb1(:,i) = slave_acc(zeros(3,1), wb1(:,i), ab1(:,i), r);
    end
  else
    zN = zeros(1,N);
    switch axs
      case 1
	accb1 =  cat(1, zN, -wb.^2*r, ab*r);
      case 2
	accb1 =  cat(1, ab*r, zN, -wb.^2*r);
      case 3
	accb1 =  cat(1, -wb.^2*r, ab*r, zN);
    end 
  end
  accb1 = accb1 + repmat(g, 1, N);
  wb1 = qrot(wb1,q12);
  accb1 = qrot(accb1,q12);

function do_unit_test()
  disp("Unit test for function track_rotation_B")
  debug = 0;

  loaddata = 0; %% Loads tracked data and generates movement from this
  
  if loaddata
    dta = load("./mocapdata/mc200_3_smooth");
    Nfr = size(dta.w,2);
    dt = 1/125;
    g = 9.82*[0;0;1];
    w0 = 10;
    r = 0.65;
    
    pL = dta.y(1:3, 1);
    pR = dta.y(4:6, 1);
    paddelc = 0.5*(pL+pR);
    vLR = pR-pL;
    vLR = vLR / norm(vLR);
    p1 = 0.65*vLR;
    p2 = -0.65*vLR;
    q12 = quaternion(randn(3,1), pi/3);
    

    wb1 = dta.w;
    accb1 = zeros(3, Nfr);
    accb2 = zeros(3, Nfr);

    for i=1:Nfr
      acenterb = qtransv(dta.a(:,i), qinv(dta.q(:,i)));
      accb1(:,i) = slave_acc(acenterb, wb1(:,i), dta.alpha(:,i), p1) + g;
      accb2(:,i) = slave_acc(acenterb, wb1(:,i), dta.alpha(:,i), p2) + g;
    end

    accb2 = qrot(accb2, q12);
    wb2 = qrot(wb1, q12);

    q0 = dta.q;
    q0inv = zeros(4,Nfr);
    for i=1:Nfr
      q0inv(:,i) = qinv(q0(:, i));
    end

    complete_revolutions = [25 176 325 469 615 763 910 1054 1201 1352 \
			    1506 1657 1803];

  else
    %% Create data of two nodes rotating about a fixed axis going through
    %% the origin and in the z-direction.
    %% Let the first imu be located at [1;0;0] initially and the second at
    %% [2;0;0]. The second is rotated about its x-axis 
    
    cycles = 6
    N = 100;
    w0 = 200 * pi / 180; % rad/s
    axes = [3,2];
    dt = 0.01;
    r = 0.6;
    r2 = -r;
    g = 9.82*randn(3,1); %% The g vector
    %%g = 9.82*zeros(3,1); %% The g vector
    d12 = 1;
    q12 = quaternion(randn(3,1), pi/3);
    %%q12 = [0 0 0 1];
    q11 = [0 0 0 1]; % no rotation
    
    p1 = [r;0.5*r;0];
    p2 = [r2;0.5*r;0];

    [wb1z, accb1z] = generate_data(N, w0, dt, 3, p1, \
				   q11, g);
    [wb1y, accb1y] = generate_data(N, w0, dt, 2, p1, \
					   q11, g);
    
    [wb2z, accb2z] = generate_data(N, w0, dt, 3, p2, \
				   q12, g);
    [wb2y, accb2y] = generate_data(N, w0, dt, 2, p2, \
				   q12, g);

    wb1 = repmat(cat(2, wb1z, wb1y), 1, cycles);
    wb2 = repmat(cat(2, wb2z, wb2y), 1, cycles);
    accb1 = repmat(cat(2, accb1z, accb1y), 1, cycles);
    accb2 = repmat(cat(2, accb2z, accb2y), 1, cycles);
    
    Nfr = size(wb1,2);
    complete_revolutions = cat(2, 1, (N-1:N:Nfr));

    %%  The orientation of the body, computed by integrating the
    %%  error-free data using the trapezoidal rule
    q0 = zeros(4,Nfr);
    q0inv = zeros(4,Nfr);
    q00 = [0;0;0;1];
    q0(:,1) = q00;
    q0inv(:,1) = q00;
    for i=2:Nfr
      q00 = qmult(q00, qexp(dt*0.5*(wb1(:,i-1)+wb1(:,i))));
      q0(:,i) = q00;
      q0inv(:,i) = qinv(q00);
    end
  end %% if loaddata

  %% Generate state sequence, and simulate sensor output
%%  xsim = cat(1, wb1, ...
%%	     (centraldiff(wb1', 1/dt))', ...
%%	     accb1, ...
%%	     zeros(3,Nfr), ...
%%	     repmat(p120, 1, Nfr));

%%  ysim = zeros(12,Nfr);
%%  for fr = 1:Nfr
%%    [rsl, ysim(:,fr)] = residual_twoimu(zeros(12,1), xsim(:,fr), q12);
%%  end
  
  cycles = length(complete_revolutions) -1;
  p120 = p2 -p1;
  q120 = q12;
  q12 = qmult(q12, quaternion(randn(3,1),pi*10/180)); % Accuracy in apriori angle measurement
  

  time1 = (1:size(wb1, 2))*dt;
  time2 = (1:size(wb2, 2))*dt;

  

  SNR = 10;
  bandwidth = 0.2;
  nreps = 1;
  fixp12 = 0;
  fixq12 = 0;
  q21 = qinv(q120);

  errangleB = zeros(length(time1), 1);
  errangleA = zeros(length(time1), 1);

  rmseA = zeros(1, length(time1));
  rmseB = zeros(1, length(time1));
  rmse_ref = zeros(1, length(time1));
  [bf,af] = butter(4, 0.25);
  meanwA = zeros(1, length(time1));
  meanwB = zeros(1, length(time1));
  meanw_ref = zeros(1, length(time1));
  errqA = zeros(1, length(time1));
  errqB = zeros(1, length(time1));
  errq_ref = zeros(1, length(time1));

  errrevolB = zeros(cycles,1);
  errrevolA = zeros(cycles,1);



  for rep = 1: nreps
    p12 = p2-p1 + randn(3,1)*0.05; % Apriori estimate of pos of sensor 2 wrt sensor 1
    q12 = qmult(q120, quaternion(randn(3,1),pi*10/180)); % Accuracy in apriori angle measurement

    eacc = randn(size(accb1))'*r*w0^2/SNR;
    eacc2 = randn(size(accb2))'*r*w0^2/SNR;
    ew = randn(size(wb1))'*w0/SNR;
    ew2 = randn(size(wb2))'*w0/SNR;
    [b,a] = butter(4,0.8);
    eaccf = eacc;
    eaccf2 =eacc2;
    ewf = ew;
    ewf2 = ew2;
    %%eaccf = filtfilt(b,a,eacc);
    %%eaccf2 = filtfilt(b,a,eacc2);
    %%ewf = filtfilt(b,a,ew);
    %%ewf2 = filtfilt(b,a,ew2);

%  tic();
    if fixp12
      p12 = p120+randn(3,1)*0.005; % Assuming known to within 5mm accuracy
    end
    if fixq12
      q12 = q120;
    end

    [x, T, qB, q12B] = track_rotation_B(time1, wb1+ewf', accb1+eaccf', ...
					time2, wb2+ewf2', accb2+eaccf2',...
					bandwidth, q12, p12, fixp12, fixq12);
%  toc()
%  tic();
%  [xs, Ts, qes, q12es] = track_rotation_B(time1, ysim(1:3,:)+ewf', ...
%				     ysim(4:6,:)+eaccf', ...
%				     time2, ysim(7:9,:)+ewf2', ysim(10:12,:)+eaccf2',...
%				     bandwidth, q12, p12);
%  toc()

%  tic();
  [xA, TA, qA, q12A] = track_rotation_A(time1, wb1+ewf', ...
					time2, wb2+ewf2', ...
					bandwidth, q12, fixq12);
%  toc()

  %% Compute the RMS error in estimated angular velocity
  %% Reference: plain filtering of the data
  wb1f = filtfilt(bf,af, (wb1+ewf')')';
  rmse_repref = sqrt(mean((wb1f(1:3,:)-wb1).^2));

  rmse_repA = sqrt(mean((xA(1:3,:)-wb1).^2));
  rmse_repB = sqrt(mean((x(1:3,:)-wb1).^2));
  rmseA = rmseA + 1/rep*(rmse_repA - rmseA);
  rmseB = rmseB + 1/rep*(rmse_repB - rmseB);  
  rmse_ref = rmse_ref + 1/rep*(rmse_repref - rmse_ref);

  meanwA = meanwA + 1/rep*(mean(xA(1:3,:) - wb1) - meanwA);
  meanwB = meanwB + 1/rep*(mean(x(1:3,:) - wb1) - meanwB);
  meanw_ref = meanw_ref + 1/rep*(mean(wb1f(1:3,:) - wb1) - meanw_ref);

  q0_ref = [0;0;0;1];
  for i=1:length(T)
    eqf = qmult(q12B(:,i), q21);
    errangleB(i) = errangleB(i) + 1/rep*(2*acos(eqf(4))*180/pi - errangleB(i));
    eqf = qmult(q12A(:,i), q21);
    errangleA(i) = errangleA(i) + 1/rep*(2*acos(eqf(4))*180/pi - errangleA(i));
    eqf = qmult(qB(:,i), q0inv(:,i));
    errqB(i) = errqB(i) + 1/rep*(2*acos(eqf(4))*180/pi - errqB(i));
    eqf = qmult(qA(:,i), q0inv(:,i));
    errqA(i) = errqA(i) + 1/rep*(2*acos(eqf(4))*180/pi - errqA(i));
    if (i>1)
      q0_ref = qmult(q0_ref, qexp(dt*0.5*(wb1f(:,i-1)+wb1f(:,i))));
      %%q0_ref = qmult(q0_ref, qexp(dt*(wb1f(:,i-1))));
    end
    eqf = qmult(q0_ref, q0inv(:,i));
    errq_ref(i) = errq_ref(i) + 1/rep*(2*acos(eqf(4))*180/pi - errq_ref(i));
  end

  %% Compute the error in orientation after a full revolution (accumulated)
  q0A = qA(:,complete_revolutions(1));
  q0B = qB(:,complete_revolutions(1));
  %%keyboard
  for i = 2:length(complete_revolutions)
    endr = complete_revolutions(i);
    eqf = qmult(qB(:,endr), q0B);
    eqa = 2*acos(eqf(4))*180/pi;
    if (eqa > 180)
      eqa = eqa - 360;
    end
    errrevolB(i-1) = errrevolB(i-1) + 1/rep*(abs(eqa) - errrevolB(i-1));
    q0B = qinv(qB(:,endr));

    eqf = qmult(qA(:,endr), q0A);
    eqa = 2*acos(eqf(4))*180/pi;
    if (eqa > 180)
      eqa = eqa - 360;
    end
    errrevolA(i-1) = errrevolA(i-1) + 1/rep*(abs(eqa) - errrevolA(i-1));
    q0A = qinv(qA(:,endr));
  end

end


%  keyboard

  figure(1)
  clf
  plot(T(1:end-3),x(1:3, 1:end-3)', 'linewidth', 2);
  hold on
  plot(time1, wb1')
  plot(time1, (wb1' + ewf))
  plot(T(1:end-3),xA(1, 1:end-3)', 'c', 'linewidth', 2);
  plot(T(1:end-3),xA(2, 1:end-3)', 'm', 'linewidth', 2);
  plot(T(1:end-3),xA(3, 1:end-3)', 'y', 'linewidth', 2);

  
  figure(2)
  clf
  plot(T(1:end-3),x(7:9, 1:end-3)', 'linewidth', 2);
  hold on
  plot(time1, accb1')
  plot(time1, (accb1' + eaccf))
 
  figure(3)
  clf
  frind = (1:length(errangleB));
  plot(frind(1:end-3), errangleB(1:end-3), 'r')
  hold on
  plot(frind(1:end-3), errangleA(1:end-3), 'm')
  
  box off
  set(findobj(gca, 'type', 'line'), 'linewidth', 3)
  set(findobj(gcf, '-property', 'fontsize'), 'fontsize', 14)
  ylabel('Degrees', 'fontsize', 14)
  xlabel('Sample number', 'fontsize', 14)
  legend('Algorithm B', 'Algorithm A')
  
  respth = date;
  mkdir(respth);

  print(fullfile(respth,"q12errorABsa.tex"), "-depslatexstandalone")
  %%print("q12errorAB.pdf", "-dpdf")
  %%print("q12errorAB.eps", "-depsc2")


  figure(4)
  clf
  cind = (1:length(errrevolB));
  plot(cind, errrevolB, 'rx', 'markersize', 16)
  hold on
  plot(cind, errrevolA, 'mx', 'markersize', 16)
  
  box off
  set(findobj(gca, 'type', 'line'), 'linewidth', 3)
  set(findobj(gcf, '-property', 'fontsize'), 'fontsize', 14)
  ylabel('Degrees', 'fontsize', 14)
  xlabel('Revolution', 'fontsize', 14)
  legend('Algorithm B', 'Algorithm A', 'location', 'north')
  set(gca, 'xlim', [0.5 length(errrevolB)+0.2])
  set(gca, 'xtick', [1 2 3 4])
  
  print(fullfile(respth,"qerrorABsa.tex"), "-depslatexstandalone")
  %print("qerrorAB.eps", "-depsc2")

  figure(5)
  clf
  frind = (1:length(rmseB));
  subplot(211)
  plot(frind(1:end-3), rmseB(1:end-3)*180/pi, 'r')
  hold on
  plot(frind(1:end-3), rmseA(1:end-3)*180/pi, 'm')
  plot(frind(1:end-3), rmse_ref(1:end-3)*180/pi, 'b')
  
  box off
  set(findobj(gca, 'type', 'line'), 'linewidth', 3)
  set(findobj(gcf, '-property', 'fontsize'), 'fontsize', 14)
  ylabel('RMS error (degrees/s)', 'fontsize', 14)
  %%xlabel('Sample number', 'fontsize', 14)
  legend('Algorithm B', 'Algorithm A')

  subplot(212)
  plot(frind(1:end-3), meanwB(1:end-3)*180/pi, 'r')
  hold on
  plot(frind(1:end-3), meanwA(1:end-3)*180/pi, 'm')
  plot(frind(1:end-3), meanw_ref(1:end-3)*180/pi, 'b')
  
  box off
  set(findobj(gca, 'type', 'line'), 'linewidth', 3)
  set(findobj(gcf, '-property', 'fontsize'), 'fontsize', 14)
  ylabel('Average residual (degrees/s)', 'fontsize', 14)
  xlabel('Sample number', 'fontsize', 14)
  legend('Algorithm B', 'Algorithm A')

  print(fullfile(respth,"w_rmseAB.tex"), "-depslatexstandalone")

  figure(6)
  clf
  cind = (1:length(errqA));
  plot(cind, errqB, 'r', 'linewidth', 4)
  hold on
  plot(cind, errqA, 'm', 'linewidth', 4)
  plot(cind, errq_ref, 'b', 'linewidth', 4)
  
  box off
  set(findobj(gca, 'type', 'line'), 'linewidth', 3)
  set(findobj(gcf, '-property', 'fontsize'), 'fontsize', 14)
  ylabel('Degrees', 'fontsize', 14)
  xlabel('Sample number', 'fontsize', 14)
  legend('Algorithm B', 'Algorithm A', 'LP filter', 'location', 'north')
    
  print(fullfile(respth, "qerrABsa.tex"), "-depslatexstandalone")
  

  figure(7)
  clf
  accb2sim = zeros(3, size(x,2));
  for i=1:size(x,2)
    [res, ysim] = residual_twoimu(zeros(12,1), x(:,i), q120, p120, q12);
    accb2sim(:,i) = ysim(10:12);
  end
  cind = (1:size(accb2sim,2));
  plot(cind, accb2sim', 'linewidth', 3)
  hold on
  plot(cind, accb2')

  figure(8)
  clf
  if (size(x,1)==15)
    plot((x(13:15,:)-repmat(p120, 1, size(x,2)))')
  end

  save("-binary", fullfile(respth, "results"))

 keyboard

