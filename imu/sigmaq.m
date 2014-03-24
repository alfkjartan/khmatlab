function [X, Wm, Wc] = sigmaq(QP,M,qi)
%% X = sigmaq(QP, M, qindstart)
%% Computes sigma points for the system with quaternion representation
%% of rotation in elements qindstart:qindend of the state vector.
%%
%% Input
%%   QP    ->   The sum of the covariance matrices QP = Q+P
%%   M     ->   The current state estimate (19x1).
%%   qindstart    ->   Row index where quaternion starts. [] means no quaternion
%%
%% Output
%%   X     <-   sigma points (19 x 37)

%% Revisions
%% 2012-04-05  If qi is not given, then assume no quaternion. If qi is
%%             empty, assume quaternions in the last elements. Determine
%%             from the sizes of the arguments.
if (nargin == 0)
  do_unit_test();
else
  
  if (nargin < 3)
    quats = 0;
  else
    quats = size(M,1) - size(QP,1);
  end

  %[L,flag] = chol(QP);
  [L] = chol(QP);
  A = L';
  %[L,D] = mchol(QP);
  % A = L*sqrt(D);

  n = size(M,1);

  alpha=[];
  beta=[];
  kappa=[];

  m = size(A,1);
  mp1 = 2*m+1;

  [Wm,Wc,c] = ut_weights(m);

  X = zeros(n, 2*m + 1);
  X(:,1) = M;
  AA = sqrt(c)*[A -A];

  if (quats == 0)
    MM = repmat(M,1, 2*m);
    X(:,2:mp1) = MM + AA;
  else
    if isempty(qi)
      %% quaternions in first elements
      qinds = 1:quats*4;
      nonqinds = setdiff(1:n, qinds);
      X(nonqinds,2:mp1) = repmat(M(nonqinds), 1, 2*m) + AA(nonqinds,:);
      
      %% Taking care of the quaternion part
      for i=1:quats
	qii = (i-1)*4+1:i*4;
	wii = (i-1)*3+1:i*3;
	Mq = M(qii);
	for j=2:(size(AA,2)+1)
	  w = AA(wii,j-1);
	  X(qii,j) = qmult(Mq,qexp(w*0.5));
	end % for j
      end % for i
    else  % qi not empty
      %% Explicit index of quaternion
      qinds = qi:qi+3;
      wii = qi:qi+2;
      nonqinds = setdiff(1:n, qinds);
      nonwinds = setdiff(1:m, wii);
      X(nonqinds,2:mp1) = repmat(M(nonqinds), 1, 2*m) + AA(nonwinds,:);
      
      %% Taking care of the quaternion part
      Mq = M(qinds);
      for j=2:(size(AA,2)+1)
	w = AA(wii,j-1);
	X(qinds,j) = qmult(Mq,qexp(w*0.5));
      end % for j
    end  % if isempty(qi)
  end % if (quats == 0) 
end 

function do_unit_test()
  tol = 1e-12;

  disp('Unit test of function sigmaq')

  Sx = randn(18,18);
  Sxm = max(max(Sx(10:12,10:12)));
  Sx(10:12,10:12) = Sx(10:12,10:12) / Sxm * 5/180 * pi;
  P = eye(18)*(10/180*pi)^2 + Sx'*Sx;
  Q = eye(18)*0.01;
  PQ = P+Q;
  M = randn(19,1);
  M(10:13) = quaternion(randn(3,1),randn(1));

  [X,Wm, Wc] = sigmaq(PQ,M,10);

  [xm,S,EE,WW] = sigma_mean_q(X, 10, Wm,Wc);

  if (norm(xm-M) > tol)
    disp('Test1: Failed')
    disp('Unexpected result:')
    disp('[M xm] ='), disp(cat(2, M, xm))
  else
    disp('Test1: OK')
  end
  

  if (norm(PQ-S) > tol)
    disp('Test2 Failed')
    disp(sprintf("Unexpected result, norm %0.8f > tol", norm((P+Q)-S)))
    disp('P+Q-S = '), disp(PQ-S)
  else
    disp('Test2: OK')
  end

  %keyboard
  

  [xm,ym, Sx, Sy, Cxy] = sigma_cov_q(X, X, 10, 10, Wm,Wc);
  if (norm((P+Q)-Cxy) > tol)
    disp('Test3: Failed')
    disp(sprintf("Unexpected result, norm %0.8f > tol", norm((P+Q)-Cxy)))
    disp('P+Q = '), disp(P+Q)
    disp('Cxy = '), disp(Cxy)
  else
    disp('Test3: OK')
  end

  M = zeros(13,1);
  M(7:10) = quaternion([0;0;1],pi/4);
  %M(7:10) = quaternion([0;0;1],0);
  
  P = eye(12)*rand(1);
  Q = eye(12)*rand(1);
  %P = eye(12)*1e-6;
  %Q = eye(12)*1e-6;

  P(7:9,7:9) = eye(3)*(pi/180*7)^2;
  Q(7:9,7:9) = eye(3)*(pi/180*15)^2;

  [X, Wm, Wc] = sigmaq(P+Q,M,7);
  
  %%p00 = randn(3,4);
  p0 = cat(2, eye(3), -eye(3));
  d0 = [1;1;1];
  Y = observe_rb(X,p0,d0,7);
  p00 = p0 + repmat(d0, 1, 6);
  
  [xm,ym, Sx, Sy, Cxy] = sigma_cov_q(X, Y, 7, [], Wm,Wc);
  
  p1 = reshape(ym,3,size(p0,2));
  d1 = mean(p1,2);
  if (norm(d1-d0) > tol)
    disp('Test4: Failed')
    disp(sprintf("Unexpected result, norm %0.8f > tol", norm(d1-d0)))
    disp('Expected '), disp(d0)
    disp('Found '), disp(d1)
  else
    disp('Test4: OK')
  end

  
  
 figure(3)
  clf
  plot3(p1(1,:), p1(2,:), p1(3,:), 'b*', 'markersize', 14)
  hold on
  plot3(0,0,0,'r*', 'markersize', 16)
  plot3(p00(1,:), p00(2,:), p00(3,:), 'g*', 'markersize', 14)
  plot3(Y(1,:), Y(2,:), Y(3,:), 'm*', 'markersize', 12)
  plot3(d0(1), d0(2), d0(3), 'c*', 'markersize', 12)
  %set(gca,'ylim',[-3 3])
  %set(gca,'xlim',[-3 3])
  %set(gca,'zlim',[-3 3])
  xlabel('x', 'fontsize', 14)
  ylabel('y', 'fontsize', 14)
  zlabel('z', 'fontsize', 14)
  %keyboard

  P = eye(6)*rand(1);
  P(4:6, 4:6) = eye(3)*(10/180*pi)^2;

  w = randn(3,1);
  w = w / norm(w);
  th = pi/3; % 60 deg
  M = cat(1, 1,2,3,(quaternion(w,th))');

  qi = 4;
 
  X = sigmaq(P,M,qi);

  Ws = zeros(3,7);
  for i=1:7
     [v,th] = quaternion(X(qi:qi+3,i));
      Ws(:,i) = v*th;
  end

 figure(2)
clf
plot3(Ws(1,:), Ws(2,:), Ws(3,:), '*', 'markersize', 14);
hold on
plot3([0 Ws(1,1)], [0 Ws(2,1)], [0 Ws(3,1)], 'r', 'linewidth', 2)
for i=1:7
   plot3([0 Ws(1,i)], [0 Ws(2,i)], [0 Ws(3,i)])
	       end
plot3([0 Ws(1,1)], [0 Ws(2,1)], [0 Ws(3,1)], 'r', 'linewidth', 3)
axis equal


  Sx = randn(18,18);
  Sxm = max(max(Sx(10:12,10:12)));
  Sx(10:12,10:12) = Sx(10:12,10:12) / Sxm * 10/180 * pi;
  %P = eye(18)*(10/180*pi)^2 + Sx'*Sx;
  P = eye(18)*(10/180*pi)^2;
  Q = eye(18)*0.01;
  PQ = P+Q;
  M = randn(19,1);
  M(10:13) = quaternion(randn(3,1),randn(1));

  [X,Wm, Wc] = sigmaq(PQ,M,10);

  [Z] = observe_imu(X);
  
  
  XX = dynamics_rbb(X,0.01);

  [Xm, Zm, Px, Pz, Pxz] = sigma_cov_q(X, Z, 10, [], Wm, Wc);

  Pv = Pz + eye(size(Pz,1));
  K = Pxz / Pv;


  keyboard

 


