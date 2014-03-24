%% Script generate_test_asc

%% Kjartan Halvorsen
%% 2013-01-29

fname = 'test.asc';
fid = fopen(fname, 'w')

%% Write a simple header
fprintf(fid, '\tfile name: test.asc\tdate/time 29.01.13 11.30\n');
fprintf(fid, 'sensor type: S2073_11-ROUND-200KPa\n');
fprintf(fid, 'total time [secs]: 12.840\ttime per frame [secs]: 0.004\tscanning rate [Hz]:  250\n');
fprintf(fid, 'pressure values in kPa\n\n');
fprintf(fid, 'electrical:\n');
fprintf(fid, 'time[secs]\n');

%% Now write the data
dta = zeros(3210, 65);
time = linspace(0, 12.840-0.004, 3210);
dta(:,1) = time;

dta(1:800, 2:4:end) = 10;
dta(801:1600, 3:4:end) = 20;
dta(1601:2400, 4:4:end) = 30;
dta(2401:3210, 5:4:end) = 40;

dlmwrite(fid, dta, '-append', 'delimiter', '\t');

fclose(fid);


