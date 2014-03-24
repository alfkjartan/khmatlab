%% Script for validation of sensor node imu against motion capture.
%% Data collected on 2012-04-26

%% Kjartan Halvorsen
%% 2012-05-04

synchronize = 1;
checkAcc = 0;

load d150w120bpm 
%% loads dennis150w120bpm, a struct containing prn and pln
%% prn and pln in turn contains the fields
%%   raw
%%   time
%%   forceBin
%%   accBin
%%   gyroBin
%%   stopwatch
pd = dennis150w120bpm;


%% Load the mocap data
dtapth = \
    '/home/kjartan/Dropbox/IJSE-Cal/data/sth-kajak-mocap/Data/dennis';
dtafile = fullfile(dtapth, '150W120bpm.tsv');


md = openmocapfile('', dtafile);
sfreq =  str2double(getvalue(md{1}, 'FREQUENCY'));


if (synchronize)
  %% time offset of same event. Found by visual inspection
  clap1 = extractmarkers(md, {'clap1'});
  clap2 = extractmarkers(md, {'clap2'});
  clap1(find(clap1==0)) = NaN;
  clap2(find(clap2==0)) = NaN;

  v = clap2 - clap1;
  vclap2 = centraldiff(clap2, 1);
  vv = sum(v.*vclap2,2);

  figure(1)
  clf
  plot(vv)
  accL = pd.pln.accBin;
  accR = pd.prn.accBin;

  accL = sum(accL.^2, 2);
  accR = sum(accR.^2, 2);

  figure(2)
  clf
  plot(accL)


  figure(3)
  clf
  plot(accR)

end

%% mocap offset---
fr0_mc = 2074;
%%-----------------

%% PLN offset---
fr0_pln = 7635;
%%-----------------

%% PRN offset---
fr0_prn = 8106;
%%-----------------


%% sync between nodes:
disp('Time pln at clap:')
disp(pd.pln.time(fr0_pln))

disp('Time prn at clap:')
disp(pd.prn.time(fr0_prn))

%if synchronize
  %% Look at the data
  nfr = 500;

  figure(3)
  clf
  t = (-50:500)*0.01;
  plot(t, accL(fr0_pln-50:fr0_pln+nfr), 'r')
  hold on
  plot(t, accR(fr0_prn-50:fr0_prn+nfr), 'g')
  
  figure(4)
  clf
  plot(t/2, clap2(fr0_mc-50:fr0_mc+nfr,3))

%end

%% Calibration data. Linear model for all sensors:
%% SI_unit = cal*(binary-offset);

%% Freescale accel: 
mV_per_g = 180;
mV_per_bit = 3.3*1000/1024;
%%g_per_bit = 1/mV_per_g*mV_per_bit;
g_per_bit = 1/mV_per_g*mV_per_bit*0.6; % Including fudge factor
apb = 9.825*g_per_bit;

%% Invensense IDG-500 gyro
mV_per_deg_per_s = 2.0*3.3/3.0;
%mV_per_bit = 3.3*1000/1024;
rad_per_s_per_bit = (pi/180.0) * 1.0/mV_per_deg_per_s * mV_per_bit;
rpspb = rad_per_s_per_bit;

gyro0_L = [420 423 417]; % bit
gyroCal_L = rpspb*eye(3); %% rad/s/bit

gyro0_R = [434 421 413]; % bit
gyroCal_R = rpspb*eye(3); %% rad/s/bit

acc0_L = [512 512 512];
accCal_L = diag([apb apb apb]); %% m/s/s/bin

acc0_R = [512 512 512];
accCal_R = diag([apb apb apb]); %% m/s/s/bin


%% Apply calibration and resample (interpolate) sensor data
timeCalL = ( pd.pln.time(fr0_pln:end) - pd.pln.time(fr0_pln) ) *1.024;
nfrLc = timeCalL(end)*100;
nfrL = length(timeCalL);
pd.pln.timeCal = linspace(timeCalL(1), timeCalL(end), nfrLc+1);

pd.pln.accCal = interp1(timeCalL, ( pd.pln.accBin(fr0_pln:end,:) - \
		 repmat(acc0_L, nfrL, 1) ) * accCal_L, pd.pln.timeCal);

pd.pln.gyroCal = interp1(timeCalL, ( pd.pln.gyroBin(fr0_pln:end,:) - \
		 repmat(gyro0_L, nfrL, 1) ) * gyroCal_L, pd.pln.timeCal);

timeCalR = ( pd.prn.time(fr0_prn:end) - pd.prn.time(fr0_prn) ) *1.024;
nfrRc = timeCalR(end)*100;
nfrR = length(timeCalR);
pd.prn.timeCal = linspace(timeCalR(1), timeCalR(end), nfrRc+1);
pd.prn.accCal = interp1( timeCalR, ( pd.prn.accBin(fr0_prn:end,:) - \
		 repmat(acc0_L, nfrR, 1) ) * accCal_R, pd.prn.timeCal);
pd.prn.gyroCal = interp1( timeCalR, ( pd.prn.gyroBin(fr0_prn:end,:) - \
		 repmat(gyro0_R, nfrR, 1) ) * gyroCal_R, pd.prn.timeCal);


%% Filter, then resample marker data
[b,a] = butter(4, 40/100);
md{2} = filtfilt(b,a,md{2});
nfr_mc = size(md{2},1) - fr0_mc + 1;
md{2} = resample(md{2}(fr0_mc:end,:), 1,2);
time_mc = (0:size(md{2},1)-1)*0.01;

%% Determine start of cycles
%% Take this to be peaks in paddel markers x-position since force data are missing.
paddelR = extractmarkers(md, {'paddel-R'});
paddelL = extractmarkers(md, {'paddel-L'});

cf_lp = 10/50; % Cutoff freq of lp filter
thr = 0.35; % Threshold for peak values
peak_width = 30; % frames
peaksR = findpeaks(paddelR(:,1), cf_lp, thr, peak_width);
peaksL = findpeaks(paddelL(:,1), cf_lp, thr, peak_width);

%% Compute the movement of the paddle.
%% The results is expressed in a coordinate system which is local
%% to each of the nodes. Reference position is coordinate positions at
%% start of each cycle

tracking_markers = {'paddel-L'
		    'paddel-R'
		    'lateral-R'
		    'front-R'
		    'back-R'
		    'lateral-L'
		    'front-L'
		    'back-L'};

gg = getMotion(md, tracking_markers);

bf_fcn = 'sensor_frame'; % Function that is called to generate local
			 % coord sys
bf_markersR.origin = {'back-R', 'front-R', 'lateral-R'}; % Centroid is used
bf_markersR.posX = 'paddel-R'; 
bf_markersR.negX = 'paddel-L'; 
bf_markersR.posY = 'front-R'; 
bf_markersR.negY = 'back-R'; 
bf_markersR.blade = 'paddel-R'; 

[gR, bladeR] = paddel_motion(gg, md, peaksR, bf_fcn, bf_markersR);
bladeposR = mean(bladeR(peaksR,:))'
[gRimu, bladeRimu] = paddel_motion_imu(pd.prn.accCal, pd.prn.gyroCal, pd.pln.accCal, pd.pln.gyroCal, peaksR, bladeposR);

keyboard

bf_markersL.origin = {'back-L', 'front-L', 'lateral-L'}; % Centroid is used
bf_markersL.posX = 'paddel-L'; 
bf_markersL.negX = 'paddel-L'; 
bf_markersL.posY = 'front-L'; 
bf_markersL.negY = 'back-L'; 
bf_markersL.blade = 'paddel-L'; 

[gL, bladeL] = paddel_motion(gg, md, peaksL, bf_fcn, bf_markersL);
bladeposL = mean(bladeL(peaksL,:))'
[gLimu, bladeLimu] = paddel_motion_imu(pd.pln.accCal, pd.pln.gyroCal, \
				       pd.prn.accCal, pd.prn.gyroCal, \ 
				       peaksL, bladeposL);


figure(7)
clf
plot(paddelR(:,1))
hold on
yl = get(gca,'ylim');
for p=peaksR
  plot([p p], yl, 'r')
end

%keyboard

figure(5)
clf
plot(pd.prn.timeCal, pd.prn.accCal)

figure(6)
clf
plot(pd.pln.timeCal, pd.pln.accCal)



%% Test orientation of right node.
frontR = extractmarkers(md, {'back-R'}); %% Seems to be swopped in this
%% file
backR = extractmarkers(md, {'front-R'});
frontL = extractmarkers(md, {'front-L'});
backL = extractmarkers(md, {'back-L'});

vyR = frontR-backR;
vyL = frontL-backL;

vynrm = mean(sqrt(sum(vyR.^2, 2)));
figure(7)
clf
plot(vyR/vynrm)

vynrm = mean(sqrt(sum(vyL.^2, 2)));
figure(8)
clf
plot(vyL/vynrm)
 
%% Check acceleration data

%% Rotation matrix taking coordinates in right body frame to spatial
%% frame in first frame (reference pos)

ex = paddelR(1,:)-paddelL(1,:);
ex = ex' / norm(ex);
ey = backR(1,:) - frontR(1,:);
ey = ey' - (ey*ex)*ex;
ey = ey / norm(ey);
ez = cross(ex,ey);

gb = cat(2, ex, ey, ez);

if (checkAcc) 
  gR = getMotion(md, {'front-R'; 'back-R'; 'lateral-R'; 'paddel-R'; 'front-L'; 'back-L'; 'lateral-L'; 'paddel-L'});

  pd.prn.accS = rotate(pd.prn.accCal, gR, Rgb);

  %% Compute the acceleration of the sensornode centroid
  lateralR = extractmarkers(md, {'lateral-R'});
  centroidR = (lateralR + backR + frontR) / 3.0;
  cvelR = centraldiff(centroidR, 100);
  caccR = centraldiff(cvelR, 100);

  figure(8)
  clf
  plot(pd.prn.timeCal, pd.prn.accS(:,3), 'r')
  hold on
  plot(time_mc, caccR(:,3)-9.825, 'b')

end
figure(9)
