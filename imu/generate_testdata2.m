function [q, w, alpha, d, v, a] = generate_testdata2()
%%   [q, w, alpha, d, v, a] = generate_testdata2()
%% Generates test data based on real data from ergometer test.
%% Corresponds to a rigid body with origin in the middle of the paddle
%%
%% Output
%%   q       ->  Orientation of the paddle
%%   w       ->  Angular velocity of the paddle
%%   alpha   ->  Angular acceleration of the paddle
%%   d       ->  Position of the center of the paddle
%%   v       ->  velocity of the center of the paddle
%%   v       ->  acceleration of the center of the paddle

%% Kjartan Halvorsen
%% 2012-05-29

debug = 1;
loadexisting = 0;

dataset = 1;

switch dataset
  case 1
    dtaf = '/home/kjartan/Dropbox/khmatlab/imu/200_3.tsv';
    markers = {'L_shaft', 'R_shaft', 'L_sensor_front', 'L_sensor_back','R_shaft', \
		 'R_sensor_front', 'R_sensor_back'};
    startfr = 1;
    nfrs = 0;
  case 2
    dtaf = '/home/kjartan/Dropbox/khmatlab/imu/Dennis200W120bpm.tsv';
    markers = {'paddel-L', 'front-L', 'back-L', 'lateral-L', ...
	       'paddel-R', 'front-R', 'back-R', 'lateral-R'};
    nfrs = 2000;
    startfr = 8000;
end 

md = openmocapfile('', dtaf);
%%keyboard
bw = 1; % Bandwidth parameter
%% Track quaternion, ang vel, ang acc, displacement, vel and acc from
%% marker trajectories 

if (startfr>1)
  md{2} = md{2}(startfr:end,:);
end
if (nfrs > 0)
  md{2} = md{2}(1:nfrs,:);
end

nfrs = size(md{2},1);

%% Set missing values to nan
md{2}(find(md{2} == 0)) = nan;

[q, w, alpha, d, v, a, y, ypred, r] = track_rigidbody(md, markers, bw);

respth = "mocapdata";
mkdir(respth)

save("-binary", fullfile(respth, "mc200_3_smooth"), "q", "w", "alpha", "d", \
     "v", "a", "y", "ypred", "r")

if debug
  %% Plot residuals
  frind = (1:nfrs)';
  figure(1)
  clf
  subplot(311)
  plot(frind, y(1,:)', 'b')
  hold on
  plot(frind, ypred(1,:)', 'r')
  subplot(312)
  plot(frind, y(2,:)', 'b')
  hold on
  plot(frind, ypred(2,:)', 'r')
  subplot(313)
  plot(frind, y(3,:)', 'b')
  hold on
  plot(frind, ypred(3,:)', 'r')
  keyboard
end


