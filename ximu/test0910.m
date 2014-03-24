%% Script used to analyze test data from 2012-09-10.
%% X-imu placed on medial part of right ankle.
%% x-axis pointing down, y-axis pointing backwards and z-axis poining
%% medially.
%% The subject (KH) attached the node, then stood straight for 10
%% seconds, then walked to starbucks, stood in line, waited for coffee,
%% walked back, stood straight, took out node, stopped logging.

%% Kjartan Halvorsen
%% 2012-09-10

dtapath = "/home/kjartan/Dropbox/projekt/nvg/data/test0910/6522";
dtaset = "06522";

nlines = get_number_of_lines(fullfile(dtapath, \
				      [dtaset, "_DateTime.csv"]))

tic()
dttme1 = read_csv(fullfile(dtapath, \
				      [dtaset, "_DateTime.csv"]));
toc()

tic()
dttme = dlmread(fullfile(dtapath, [dtaset, "_DateTime.csv"]),\
		',', 1, 0);
toc()
tic()
%caldta = dlmread(fullfile(dtapath, [dtaset, "_CalInertialAndMag.csv"]),\
%		',', 1, 0);
[caldta, calheader] = read_csv(fullfile(dtapath, [dtaset, "_CalInertialAndMag.csv"]));
toc()


% The quaternion data
tic()
[qdta, qheader] = read_csv(fullfile(dtapath, [dtaset, "_Quaternion.csv"]));
toc()

