function sd = opensensordata(fname)
%%  sd = opensensordata(fname)
%% Reads in raw sensor data in fname. The values are in bits up to 1024
%% and separated by commas. The order of the columns are
%%   Time, 

%% Kjartan Halvorsen
%% 2012-05-18

dta = dlmread(fname, ",");

sd.raw = dta;
sd.time = dta(:,1);
sd.forceBin = dta(:,2);
sd.accBin = dta(:,3:5);
sd.gyroBin = dta(:,6:8);
