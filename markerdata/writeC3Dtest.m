% Script for testing the writeC3D.m function
%

% Kjartan Halvorsen
% 2004-08-19

olddir = pwd;
[pth,name] = fileparts(mfilename('fullpath'));
cd(pth);

filename_in = '..\kinematik\jc\4ms0001.c3d';
filename_out = '..\kinematik\jc\4ms0001_out.c3d';

[Markers,VideoFrameRate,AnalogSignals,AnalogFrameRate,Event, ...
 ParameterGroup,CameraInfo,ResidualError] = readC3D(filename_in);

writeC3D(Markers,VideoFrameRate, ...
	 AnalogSignals, AnalogFrameRate, ...
	 Event,ParameterGroup,CameraInfo, ...
	 ResidualError, filename_out)


%eval(['!diff --binary ', filename_in, ' ', filename_out])


[Markers2,VideoFrameRate2,AnalogSignals2,AnalogFrameRate2,Event2, ...
 ParameterGroup2,CameraInfo2,ResidualError2] = ...
    readC3D(filename_out);

cd(olddir);


