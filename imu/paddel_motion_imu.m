function [gg, pb, vb] = paddel_motion_imu(acc1, w1, acc2, w2, pks, blade)
%%  [gg, pb] = paddel_motion_imu(acc1, w1, acc2, w2, pks, blade)
%% Computes the movement of the paddel from imu data. 
%%
%% Input
%%   acc1             ->  accelerometer data from same side node
%%   w1               ->  gyro data from same side node
%%   acc2             ->  accelerometer data from contralateral node
%%   w2               ->  gyro data from contralateral node
%%   pks              ->  list of indices where cycle starts
%%   blade            ->  coordinates of blade midpoint (end of shaft)

%% Kjartan Halvorsen
%% 2012-05-11

nfrs = size(acc1,1);
gg = zeros(4,4,nfrs);
pb = zeros(nfrs, 3);
vb = zeros(nfrs, 3);

qf = [0;0;0;1]; %% Initial guess is no rotation between nodes

dt = 0.01; %% Sampling interval

for i=1:length(pks)
  p = pks(i);

  if (i < length(pks))
    stopfr = pks(i+1)-1;
  else
    stopfr = size(acc1,1);
  end
  
  bwidth = 1;
  ekf = 1;


  [d, vel, accs, q, ws, qfhat] = track_twonode(acc1(p:stopfr,:)', \
						w1(p:stopfr,:)', \
						acc2(p:stopfr,:)', \
						w2(p:stopfr,:)', \
						dt, bwidth, qf,ekf);
  qf = qfhat(:,end);

  g = mean(accs,2);
  
  ds = zeros(3,1);
  vs = -cross(ws(:,1), blade); %% Make velocity of blade point zero at
  %% start of cycle.
  for fr=p:stopfr
    ds = ds + dt*vs;
    vs = vs + dt*(accs(:,fr-p+1)-g);
    pb(fr,:) = ( ds + qtransv(blade, q(:,fr-p+1)) )';
    gg(:,:,fr) = eye(4);
    gg(1:3,4,fr) = ds;
    gg(1:3,1:3,fr) = quaternion2rotation(q(:,fr-p+1));
  end
end

