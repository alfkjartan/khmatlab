function [f_ex] = get_external_force(md, uppoint, downpoint, mass)
  %% [f_ex] = get_external_force(md, uppoint, downpoint)
  %% Computes the external force from the accelerating mass whose movement is given by 
  %% the movement along the line through the two markers uppoint and downpoint.
  %%
  %% Input
  %%   md          ->  Marker data
  %%   uppoint     ->  Name of point 1
  %%   downoint    ->  Name of point 2
  %%   mass        ->  Mass of moving body
  %% Output
  %%   f_ex        <-  External force

  %% Kjartan Halvorsen
  %% 2013-01-28

  %% Solution:
  %% 1) Marker data is assumed to be filtered.
  %% 2) Gravitity is assumed to have direction [0;0;-1].
  %% 3) Acceleration is calculated by twice differentiating the positions 

if nargin == 0
  do_unit_test()
else
  point1 = extractmarkers(md, uppoint);
  point2 = extractmarkers(md, downpoint);
  
  freq = str2double(getvalue(md{1}, 'FREQUENCY'));
  d = 0.5*(point2+point1);
  v = centraldiff(d, freq);
  a = centraldiff(v, freq);
  Nfr = size(a,1);

  g = [0, 0, -9.82];
  f_ex = repmat(g*mass, Nfr, 1) - a*mass;
end

function do_unit_test()
  pkg load signal


  datapth = '/home/kjartan/Dropbox/projekt/golf/data/validering120207';
  md = openmocapfile('', fullfile(datapth, 'weight1000_1.tsv'));

  freq = str2double(getvalue(md{1}, 'FREQUENCY'));
  [bf,af] = butter(4,2*10/freq); %% Filter parameters for low pass filter to compare with

  md{2} = filtfilt(bf,af, md{2});
  upp = 'weight_top';
  downp = 'weight_low';
  
  mass = 1;

  f_ex = get_external_force(md, upp, downp, mass);

  figure(1)
  clf
  plot(f_ex);
  title('External force')

  figure(2)
  clf
  plot(extractmarkers(md, upp), 'linewidth', 2)
  plot(extractmarkers(md, downp))
  title('Marker position')


