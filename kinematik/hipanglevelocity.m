function [left_ang, left_vel, right_ang, right_vel] = hipanglevelocity
%  [ang, vel] = hipanglevel
% Opens a tsv file, calculates the hip angle (deg) and angle
% velocity (deg/s).
% Markers used
% great troch       
% lateral knee
% shoulder
% ASIS
%

% Kjartan Halvorsen
% 2004-01-21

% Cut-off frequency for low-pass filtering of angle velocity.
cutoff = 10; %Hz

% Load a tsv file
[mdata, adata, filename] = opentsv;

vertforcename = 'Fz';

shld_l = extractmarkers(mdata, 'left_shoulder');
shld_r = extractmarkers(mdata, 'right_shoulder');
gt_l = extractmarkers(mdata, 'left_grtroch');
gt_r = extractmarkers(mdata, 'right_grtroch');
knee_l = extractmarkers(mdata, 'left_knee');
knee_r = extractmarkers(mdata, 'right_knee');
asis_l = extractmarkers(mdata, 'left_ASIS');
asis_r = extractmarkers(mdata, 'right_ASIS');

v_thigh_l = knee_l - gt_l;
v_thigh_r = knee_r - gt_r;
v_rl = asis_l - asis_r;

nfr = size(gt_l, 1);

left_ang = zeros(nfr,1);
left_vel = zeros(nfr,1);
right_ang = zeros(nfr,1);
right_vel = zeros(nfr,1);

for fr = 1:nfr
  % Left hip
  hipjc = gt_l(fr,:) - 20*v_thigh_l(fr,:)/norm(v_thigh_l(fr,:));
  e_th = knee_l(fr,:) - hipjc;
  e_th = e_th / norm(e_th);
  
  e_ub = shld_l(fr,:) - hipjc;
  e_ub = e_ub / norm(e_ub);
 
  left_ang(fr) = acos(e_th*e_ub');
  
  %  crossprod = cross(e_ub, e_th);
%  left_ang(fr) = sign(crossprod*v_rl(fr,:)')*asin(norm(crossprod));

  % Right hip
  hipjc = gt_r(fr,:) - 20*v_thigh_r(fr,:)/norm(v_thigh_r(fr,:));
  e_th = knee_r(fr,:) - hipjc;
  e_th = e_th / norm(e_th);
  
  e_ub = shld_r(fr,:) - hipjc;
  e_ub = e_ub / norm(e_ub);
  
  right_ang(fr) = acos(e_th*e_ub');

  %crossprod = cross(e_ub, e_th);
  %right_ang(fr) = sign(crossprod*v_rl(fr,:)')*asin(norm(crossprod));
end

left_vel = diff(left_ang);
right_vel = diff(right_ang);

% Filter data

samplefreq = str2double(getvalue(mdata{1}, 'FREQUENCY'));

if  ~isnan(samplefreq)
  order = 4;
  [b,a] = butter(order, min(1, 2*cutoff/samplefreq));
  left_vel = filtfilt(b,a,left_vel);
  right_vel = filtfilt(b,a,right_vel);
  dt = 1/samplefreq;
else
  dt=1;
end

left_vel = left_vel / dt;
right_vel = right_vel / dt;


% Check if exists analog data
if isempty(adata)
  hasanalog = 0;
else
  hasanalog = 1;
  afreq = str2double(getvalue(adata{1}, 'FREQUENCY'));
  vertforce = extractchannels(adata, vertforcename);
  oversample = afreq / samplefreq;
  if (oversample > 1)
    vertforce = decimate(vertforce, oversample);
  end
  % Make range < 100, for plotting together with angle velocity
  vertforce = 100* vertforce / max(vertforce);
end



figure
set(gcf, 'NumberTitle','off');
set(gcf, 'Name', filename);
timevel = dt*(1:length(left_vel));
timeang = dt*(1:length(left_ang));

tmax = max(timeang);

subplot(221)
plot(timeang, left_ang*180/pi)
ylabel('Degrees')
title('Left hip')
set(gca, 'xlim', [0 tmax])
subplot(223)
if hasanalog
  plot(timevel, left_vel*180/pi, timeang, vertforce)
else
  plot(timevel, left_vel*180/pi)
end  
ylabel('Degrees/s')
xlabel('s')
set(gca, 'xlim', [0 tmax])

subplot(222)
plot(timeang, right_ang*180/pi)
ylabel('Degrees')
title('Right hip')
set(gca, 'xlim', [0 tmax])
subplot(224)
if hasanalog
  plot(timevel, right_vel*180/pi, timeang, vertforce)
else
  plot(timevel, right_vel*180/pi)
end  
ylabel('Degrees/s')
xlabel('s')
set(gca, 'xlim', [0 tmax])


khazoomsafe('on')




  