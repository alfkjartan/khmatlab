function [pos,fr] = find_max(traj, dir, startfr, stopfr)
% Returns the frame with maximum in the given direction within
% start and stop frames

% Kjartan Halvorsen 
% 2010-03-24

if nargin==0
  % unit test
  testme;
  return
end

if nargin < 5
  stopfr = size(traj,1);
end

if size(traj,2) == 1
  sig = traj;
else
  dir = dir(:);

  sig = traj*dir;
end

startfr = max(1, startfr);

if stopfr > length(sig)
  stopfr = length(sig);
end

fr = (startfr - 1) ...
     + find(sig(startfr:stopfr) == max(sig(startfr:stopfr)));

%keyboard

pos = traj(fr,:);

function testme
w = 2*pi/100;
traj = cat(2, sin(w*(1:100)'), zeros(100,2));
dir = [1; 0; 0];

expected1 = 25;
found1 = find_max(traj,dir,1, 50);

if (expected1 == found1) 
  disp('Unit test 1 passed')
else
  disp(['Unit test 1 failed. Expected ', int2str(expected1), ...
	'. Found ', int2str(found1)])
end


expected1 = 75;
found1 = find_max(-traj,dir,50,200);

if (expected1 == found1) 
  disp('Unit test 2 passed')
else
  disp(['Unit test 2 failed. Expected ', int2str(expected1), ...
	'. Found ', int2str(found1)])
end


