function [f_ex] = do_unit_test(md,tokens, event1, event2)
  %pkg load signal


%   datapth = 'C:\Users\fredrikt\Documents\MATLAB\PlianceValidation\20130410';
%   md = openmocapfile('', fullfile(datapth, '2000_100_7SG.tsv'));

%  freq = str2double(getvalue(md{1}, 'FREQUENCY'));
%  [bf,af] = butter(4,2*10/freq); %% Filter parameters for low pass filter to compare with

%  md{2} = filtfilt(bf,af, md{2});
  upp = 'weight_top';
  downp = 'weight_low';
  
  weight = 1; %str2num(tokens{1}{1})/1000;
  mass = 0.1185 + weight;

  f_ex = get_external_force(md, upp, downp, mass);

  figure(2)
  clf
  plot(extractmarkers(md, upp), 'linewidth', 2)
  plot(extractmarkers(md, downp))
  title('Marker position')
