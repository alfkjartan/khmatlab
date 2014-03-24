function ok = test_padkah2_notw
%  ok = test_padkah2_notw
% Unit test for padkah2_notw. Returns ok if all tests passed.
%  
% This function tests not only for bugs, but also the robustness of
% the algorithm.
 
% Kjartan Halvorsen
% 2004-04-27
  
tol = 1e-12;

% Create a number of test cases
wa1 = [1;0;0];
%phi = linspace(0,7*pi/16,10);
phi = linspace(0,1*pi/16,10);
wa2 = cat(1, sin(phi), cos(phi),zeros(size(phi)));
qa = zeros(3,1);

twa1 = create_twist(wa1,qa);
twa2 = zeros(4,4,10);
for i=1:10
  twa2(:,:,i) = create_twist(wa2(:,i),qa);
end

wtmp = rand(3,1); wtmp = wtmp / norm(wtmp);
twtemp = create_twist(wtmp,zeros(3,1));
gba = expr(twtemp*rand(1,1));
Rba = gba(1:3,1:3);
wb1 = Rba*wa1;
wb2 = Rba*wa2;

qb = ones(3,1);
%qb = qa;
twb1 = create_twist(wb1,qb);
twb2 = zeros(4,4,10);
for i=1:10
  twb2(:,:,i) = create_twist(wb2(:,i),qb);
end

phi = linspace(0,1*pi/16,10);

pp = cat(1,sin(phi),zeros(size(phi)), cos(phi), ones(size(phi)));
ppb = (pp + repmat(cat(1,qb,0),1,size(pp,2)));


th1 = linspace(4*pi/16,7*pi/16,5);
th2 = linspace(4*pi/16,7*pi/16,5);

% Test for bugs
k=1;
for t1=th1
  ga1 = expr(twa1*t1);
  gb1 = expr(twb1*t1);
  for t2=th2
    ok = 1;
    for ii=1:10
      ga = ga1*expr(twa2(:,:,ii)*t2);
      gb = gb1*expr(twb2(:,:,ii)*t2);
      qqa = ga * pp; 
      qqb = gb * ppb; 
      for jj=1:10
	[tha,flag,thalta] = ...
	    padkah2_notw(qa, wa1, wa2(:,ii), pp(1:3,jj), qqa(1:3,jj)); 
	[thb,flag,thaltb] = ...
	    padkah2_notw(qb, wb1, wb2(:,ii), ppb(1:3,jj), qqb(1:3,jj)); 
	oka = assert([t1; t2],tha,tol);
	okb = assert([t1; t2],thb,tol);
	if (oka & okb)
	  %disp(['Test ', int2str(k), ' ok'])
	else
	  %disp('Alternative solution:')
	  %thaa = [tha thalta]
	  %thbb = [thb thaltb]
	end
	k=k+1;
	ok = ok & oka & okb;
      end
    end
    if ok
      disp(['th1=',num2str(t1),'   th2=',num2str(t2),'.....ok'])
    else
      disp(['th1=',num2str(t1),'   th2=',num2str(t2),'.....failed'])
    end
    
  end
end

    
  
  
  
% Test for robustness
  

function equals = assert(expected, found, tolr)
  % Checks that the found equals expected. Prints a report on the
  % screen if not.
    
  if (nargin==2)
    equals = min(min(expected == found));
  else
    equals = min(min(abs(expected - found) < tolr));
  end

  if ~equals
%    disp('Assertion failed!')
%    expected = expected
%    found = found
  end
  
  
function twm = create_twist(w,q)
  % Utility function for creating twist matrices
  % Input
  %    w     ->  Unit vector in direction of twist axis  
  %    q     ->  A point on the axis
    
  v = -cross(w,q);
  twm = hat([v;w]);
  