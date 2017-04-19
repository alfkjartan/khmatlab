function [Jsb,ndofs]=link_jacobians(tws, gsli, th, thind, thproxind, gprox, twsprox)
% Jsb=link_jacobians(tws, gsli, th)
% Returns the body jacobian for each link in km. The body jacobian J_{st}^b is defined by the 
% relationship
%          V_{st}^b = J_{st}^b \dot{\theta}
% Or in words: The linear relation between joint velocities and the segment velocity in 
% body coordinates. Column i is given by the ith instanenous joint twist relative to 
% the segment frame.  
% Input
%   tws    ->  nested array of twists
%   gsli   ->  nested array of transformations taking points in the local coordinate system of
%              each link to the static frame for th=0
%   th     ->  generalized coordinates
%   thind  ->  index into generalized coordinates one before start of coordinates 
%              for current segment
%   thproxind  ->  list of indices into generalized coordinates that this segment depends upon 
%   twsprox -> twists from root up to current segment (6 x nprox) matrix
%   gprox   -> proximal rigid transformations (4 x 4 x nprox) matrix, where the 
%              ith transformation is
%                 g_i = expr(twsprox(1), thprox(1))* ... * expr(twsprox(i), thprox(i))
% Output
%   Jsb    <-  The jacobians, a (6 x nsts x nsgms) matrix

%% Kjartan Halvorsen
% 2013-05-29

if nargin == 0
   do_unit_test();
else

  %% Total number of states
  nststotal = size(th,1);

  if nargin == 3  % Initial call
    gprox = zeros(4,4,nststotal);
    twsprox = zeros(4,4,nststotal);
    thproxind = [];
    thind = 0;
    nprox = 0;
  else
    nprox = size(twsprox, 3);
  end

  %% First the "own" Jacobian, then the branches
  mygsl = gsli{1};
  mytws=tws{1};

  nn=length(mytws);
  myinds = thind+1:thind+nn;
  thproxind = cat(2, thproxind, myinds);
  myx=th(myinds);
 
  for st=myinds
    gprox(:,:,st)=expm( mytws{st-thind}*th(st) );
    twsprox(:,:,st)=mytws{st-thind};
  end

  %%gprox
  %%twsprox

  Jsb = zeros(6, nststotal); 
  gg=gsli{1};
  for st=fliplr(thproxind)
    gg=gprox(:,:,st)*gg;
    Jsb(:,st) = adjoint_trf_inv(vee(twsprox(:,:,st)), gg);
  end

  % ------------------------------------------------
  %  The branches
  % ------------------------------------------------

  ndofs = nn;
  thind = thind + nn;
  try
  if (length(tws)>1) % Branches exist
    %% For debug  warning(['Found ', int2str(length(tws)-1), ' branches'])
    twsbr=tws(2:length(tws));
    gslbr=gsli(2:length(tws));
    ndist=0;
    nsts = size(th,1);
    for br=1:length(twsbr)
      [Jsbdist, ndofsbr] = link_jacobians(twsbr{br}, gslbr{br},...
			       th, thind, thproxind, ...
			       gprox, twsprox);
      ndofs = ndofs + ndofsbr;
      thind = thind + ndofsbr;
      Jsb = cat(3, Jsb, Jsbdist);
    end
  end
  catch 
    keyboard
  end

end


function do_unit_test()

l1 = 1;
l2 = 2;
m1 = 3;
m2 = 2;
m3 = 0.6;
m4 = 0.3;

th = [pi/4; pi/5; pi/6; pi/7];
th1 = th(1);
th2 = th(2);

sm = scara_robot_model(l1, l2, m1, m2, m3, m4);

[Mtrue, M1, M2, M3, M4] = scara_robot_inertia(sm, th2);

p1 = [0;l1/2;0];
p2 = [0;l1+l2/2;0];
p3 = [0; l1+l2;0];
p4 = [0;l1+l2;1];

I3 = eye(3);
Z3 = zeros(3,3);

tws = { { hat([0;0;0;0;0;1]) }, 
	{ { hat([l1;0;0;0;0;1]) }, 
	  { { hat([l1+l2;0;0;0;0;1]) }, 
	    { { hat([0;0;1;0;0;0]) } }}}};

gsli = {  [I3 p1; 0 0 0 1] , 
	{  [I3 p2; 0 0 0 1] , 
	  {  [I3 p3; 0 0 0 1] , 
	    {  [I3 p4; 0 0 0 1]  }}}};


Jsb = link_jacobians(tws, gsli, th);

Mlinks = cat(3, M1, M2, M3, M4);

M = zeros(4,4);
for s=1:4
  Jsbs = Jsb(:,:,s);
  M = M + Jsbs'*Mlinks(:,:,s)*Jsbs;
end

tol = 1e-12;

if norm(M - Mtrue) > tol
   disp('Test 1 failed')
   disp('Expected'), Mtrue
   disp('Found'), M
   keyboard
else
   disp('Test 1 OK')
end

