function [nd, d, vel, acc, q, w, alpha] = track_rigid(mdata, dt, \
						      bandwidth, q0)
% [nd, d, vel, acc, q, w, alpha] = track_rigid(mdata, orig, R0, dt, bandwidth)
% Tracks a rigid body from marker data using a fix interval ukf
% smoother.
% Uses the implementation of ukf in the EKF/UKF toolbox from Aalto university.
%
% Input
%    mdata      ->   marker data (nfr x 3*nmarkers).
%    dt         ->   sampling period
%    bandwidth  ->   Filter tuning parameter
% Output
%    nd         <-   a matrix with smoothed marker data (fixed
%                       interval smoother, based on all data). Use this!
%    d          <-   the estimated translation of the rigid body. The
%                    local coordinate of the RB coincides with the
%                    static lab frame for the first frame of data.
%    vel        <-   velocity of the RB, expressed in the static frame.
%    acc        <-   acceleration of the RB, expressed in the static frame.
%    q          <-   quaternion representing the rotation of the RB.
%                    Rotates a point in the first frame to subsequent positions.
%    w          <-   Angular velocity of RB. Expressed in the static
%                    frame.
%    alpha      <-   Angular accelerastion of RB. Expressed in the
%                   static frame.

% Kjartan Halvorsen
% 2012-03-20. 
%

if (nargin == 0)
  do_unit_test();
else

  if (nargin < 4)
    q0 = [0;0;0;1];
  end

  n = 19; % Or 19, the size of the state vector
  %n = 13; % Or 19, the size of the state vector

  %% The state vector is x=[d,v,acc,q,w,alpha]
  Y = mdata'; % point observations in columns
  nfr=size(Y,2);
  nmarkers = size(Y,1)/3;

				% The 6dof rigid body model
  p00 = reshape(Y(:,1), 3, nmarkers);
  d0 = mean(p00,2); % The centroid
  p0 = p00 - repmat(d0, 1, nmarkers);
  p0q = p0;
  for i=1:nmarkers
    p0q(:,i) = qtransv(p0(:,i), qinv(q0));
  end
  p0 = p0q;

  sd = 0.01; % in m standard deviation on marker positions
  %sd = 1;
  R = eye(3*nmarkers) * sd^2 / bandwidth;

				% Spectral power density of the white noise.
  %%   displ        vel       acc      quat                  w      alpha 
  if (n == 19)
    q = cat(2, [0.05 0.05 0.05],[1 1 1],[10 10 10],...
	    [1 1 1]*pi/180*10, [1 1 1], [10 10 10]).^2;
    qindstart = 10;
  else
    q = cat(2, [0.05 0.05 0.05], [1 1 1],  [1 1 1]*(pi/180*10), [1 1 1]).^2;
    qindstart = 7;
  end

  Q = diag(q);
  
  nframes = size(Y, 2);

  x_pred = zeros(n,nframes);
  x_upd = zeros(n,nframes);
  Px = zeros(n-1,n-1,nframes);

  P = Q*4;
  x = zeros(n, 1);
  
  %x(qindstart+3) = 1;
  x(qindstart:qindstart+3) = q0;
  %x(qindstart:qindstart+3) = quaternion(randn(3,1), pi/2)';

  %vs = zeros(nmarkers*3, nframes);
  %Ks = zeros(n-1, nmarkers*3, nframes);
  %Pxzs = zeros(n-1,nmarkers*3, nframes);
  try
    for fr = 1:nframes
      [x, P] = ukf_predictq(x, P, Q, qindstart, dt);
      %disp(fr)
      x_pred(:,fr) = x;
      [x, P, v, K,Pxz] = ukf_updateq(x,P,Y(:,fr), R, qindstart, p0, d0);
      x_upd(:,fr) = x;
      %Px(:,:,fr) = P;
      %vs(:,fr) = v;
      %Ks(:,:,fr) = K;
      %Pxzs(:,:,fr) = Pxz;
      %disp(fr)
    end
  catch 
    keyboard
  end

				%[xs,Ps]=utf_smooth(xM, Px, Y, f_func, Q, dt, h_func, R, p0 );

  %nd = (observe_rb(x_upd, p0, d0, qindstart))';
  nd = (observe_rbb(x_upd, p0, d0, qindstart))';

  d = x_upd(1:3,:);
  vel = x_upd(4:6,:);

  if (n==19)
    acc = x_upd(7:9,:);
    q = x_upd(10:13,:);
    w = x_upd(14:16,:);
    alpha = x_upd(17:19);
  else
    acc = [];
    q = x_upd(7:10,:);
    w = x_upd(11:13,:);
    alpha=[];
  end
end

function do_unit_test()
  disp("Unit test for function track_rigid")

  tol = 1e-10;

  N = 200;
  m = 6;
  md = repmat(randn(1,m*3)*0.1, N, 1);
  md(4,:) = md(4,:) + repmat(randn(1,3)*0.2, 1, m);
  dt = 0.01;

  w0 = randn(3,1);
  th0 = pi/6;
  %q0 = quaternion(w0,th0)';
  q0 = [0;0;0;1];

  % Ramp
  w = randn(1,3)*0.5;
  d = linspace(0,0.2,30);
  dd = cat(2, d, fliplr(d));
  for j=21:21+59
    md(j,:) = md(j,:) + repmat(w*dd(j-20), 1, m);
  end

  % Rotate grad degrees
  grad = 270;
  periods = 2;
  w = randn(3,1);
  w = w/norm(w);
  th = sin(linspace(0, 2*pi*periods, 80))* pi/180 * grad;
  for j=80:80+79
    qw = quaternion(w,th(j-80+1));
    md20 = reshape(md(j,:),3,m);
    for i=1:m
      md20(:,i) = qtransv(md20(:,i), qw);
    end
    md(j,:) = reshape(md20, 1, 3*m);
  end

  ee = randn(size(md))*0.005;
  [b,a] = butter(4,0.8);
  ef = filtfilt(b,a,ee);
  [nd, d, v, a, q, w, aa] = track_rigid(md + ef, dt, 1, q0);

  figure(1)
  clf

  fr = 1:size(nd,1);
  plot(fr, nd(:,1:3))
  hold on
  plot(fr, md(:,1), 'color', [0.4 0.4 1])
  plot(fr, md(:,2), 'color', [0.4 1 0.4])
  plot(fr, md(:,3), 'color', [1 0.4 0.4])

  plot(fr, md(:,1)+ef(:,1), 'color', [0.6 0.6 1])
  plot(fr, md(:,2)+ef(:,2), 'color', [0.6 1 0.6])
  plot(fr, md(:,3)+ef(:,3), 'color', [1 0.6 0.6])

  keyboard


