%function [acc1, w1, acc2, w2, d, q, d2, q2] = generate_testdata()
%%  [acc1, w1, acc2, w2, d, q] = generate_testdata()
%% Generates test data based on real data from ergometer test.
%%
%% Output
%%   acc1    ->  Acceleration at node 1 (3 x nframes)
%%   w1      ->  Angular velocity at node 1 (3 x nframes)
%%   acc2    ->  Acceleration at node 2 (3 x nframes)
%%   w2      ->  Angular velocity at node 2 (3 x nframes)
%%   d       ->  Path of origin of node 1 (in lab frame) (3 x nframes)
%%   q       ->  Orientation of node 1 (in lab frame) (4 x nframes)
%%   d2      ->  position of node 2 in local coord sys of node 1
%%   q2      ->  orientation of node 2 in local coord sys of node 1

%% Kjartan Halvorsen
%% 2012-03-19
debug = 0;
loadexisting = 1;

dtaf = '/home/kjartan/Dropbox/khmatlab/imu/200_3.tsv';
md = openmocapfile('', dtaf);

bw = 1; % Bandwidth parameter
%% Track quaternion, ang vel, ang acc, displacement, vel and acc from
%% marker trajectories 

dt = 1/str2double(getvalue(md{1}, 'FREQUENCY'));

l_shaft = extractmarkers(md, {'L_shaft'});
l_sensor_fr = extractmarkers(md,{'L_sensor_front'});
l_sensor_back = extractmarkers(md, {'L_sensor_back'});
r_shaft = extractmarkers(md, {'R_shaft'});
r_sensor_fr = extractmarkers(md, {'R_sensor_front'});
r_sensor_back = extractmarkers(md, {'R_sensor_back'});

shaftm = cat(2, l_shaft, l_sensor_fr, l_sensor_back, r_shaft, ...
	     r_sensor_fr, r_sensor_back);

%shaftm = shaftm(1:300,:);

%% Local coordinate system fixed to the shaft. 
%% Origin is in R_sensor_back, y-axis is in direction of the shaft,
%% rigth -> left, x axis is in the direction R_sensor_back ->
%% R_sensor_front
orig = r_sensor_back(1,:)';
ey = l_shaft(1,:) - r_shaft(1,:);
ey = ey' / norm (ey);
ex = (r_sensor_fr(1,:) - r_sensor_back(1,:));
ex = ex' - (ex*ey)*ey;
ex = ex / norm(ex);
ez = cross(ex,ey);
R0 = cat(2, ex, ey, ez);
%% R0 operated on vector v=[v1 v2 v3]' in local coordinate system will
%% give vs = v1*ex + v2*ey + v3*ez, that is
%% the same vector expressed in the global coordinate system


if loadexisting
  load paddle.dat
else
  [shaftm_filtered, d, vel, acc, q, w, alpha] = ...
      track_rigid(shaftm, dt, bw);
end

nmarkers = size(shaftm,2)/3;
p00 = reshape(shaftm(1,:), 3, nmarkers);
d0 = mean(p00,2); % The centroid

%% Transform acceleration and gyro data to local coordinate system of node 1
qsb = rot2quat(R0');% Rotates  vector in  body frame to spatial frame
qbs = qinv(qsb); % Rotates  vector in spatial frame to body frame

qq10 = qrot(q,qbs);


%% Compute acc, w
acc1 = qrot(acc,qbs);
w1 = qrot(w, qbs);

%% Track from imu data

[d1, vel1, q1] = track_imu(acc1, w1, dt, bw);

%% Determine d2, q2
%v = [0;1;0]
%th = pi/3; % 60 degrees rotation between the two blades
%q2 = quaternion(v,th);
%d2 = [0;0.5; -0.03] + qrot(q2, [0;0;0.03]);


%acc2s = acc 
 
%acc2 = qrotconj(qprod(q,q2), acc2s);
%w2 = qrotconj(qprod(q,q2), w);

if debug
  %% Plot residuals
  plotmarkers(shaftm, shaftm_filtered)

  figure(7)
  clf
  subplot(311)
  plot(w(1,:)*180/pi)
  ylabel('X Deg/s')
  subplot(312)
  plot(w(2,:)*180/pi)
  ylabel('Y Deg/s')
  subplot(313)
  plot(w(3,:)*180/pi)
  ylabel('Z Deg/s')

  figure(8)
  clf
  subplot(311)
  plot(acc(1,:)/9.82)
  ylabel('X g')
  subplot(312)
  plot(acc(2,:)/9.82)
  ylabel('Y g')
  subplot(313)
  plot(acc(3,:)/9.82)
  ylabel('Z g')

  figure(9)
  clf
  subplot(311)
  plot(d(1,:))
  ylabel('X displ m')
  subplot(312)
  plot(d(2,:))
  ylabel('Y displ m')
  subplot(313)
  plot(d(3,:))
  ylabel('Z displ m')

end





