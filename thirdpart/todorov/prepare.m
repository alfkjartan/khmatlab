%[map, info, S, R, V] = prepare(segment,mode, [dt,datS,datR,datV,ratio])
%
% Prepare model for estimation
%
% segment: kinematic model (see residual.m for details)
% mode: integer specifying which joints should be estimated
%       1- variable only; 2- constant only; 3- variable and constant
%
% dt: sampling time step
% datS: 2x1 vector specifying the positional standard deviations for
%        [variable, constant] components of the initial state
% datR: positional stdev for dynamics noise, same format as datS
% datV: scalar specifying positional stdev of sensor noise
% ratio: positional/angular stdev; note that angles are always in rad
%         while positions are in arbitrary units
%
% map: structure describing the model (see residual.m for details)
% info: structure with fields 'type' and 'spatial'; each field is a
%        vector of flags with the dimensionality of the state vector
%       type: 1- component of variable joint; 2- constant joint
%       spatial: 1- component of translation joint; 2- rotation joint
%
% S,R,V: covariance matrices computed from the optional input arguments
%         dt,datS,datR,datV,ratio
%        if S,R,V are not assigned the optional inputs can be omitted

% Copyright (C) Emanuel Todorov, 2006-2007

function [map, info, S, R, V] = ...
   prepare(segment, mode, dt, datS, datR, datV, ratio)

tpS = 1;                                     % type Slide
tpT = 2;                                     % type Translation
tpH = 3;                                     % type Hinge
tpR = 4;                                     % type Rotation
szj = [1 3 1 3];                             % size of joint/residual data

nS = size(segment,1);                        % get number of segments


%---------- construct joint map ----------------------------------------
z = zeros(nS,1);
map = struct('active',false(nS,1),...        % initialize map structure
   'aR',z,'aY',z,'aJ',z,'aP',z,'aW',z,...
   'nR',0,'nY',0,'nJ',0,'nP',0,'nW',0);
   
map.active = (bitand(segment(:,5),mode)>0);  % set active joints

sz = (segment(:,3)==tpH) + 3*(segment(:,3)==tpT) + 3*(segment(:,3)==tpR);
[map.aR, map.nR] = address(sz);              % addresses in r

sz = (segment(:,3)==tpH) + 3*(segment(:,3)==tpT) + 4*(segment(:,3)==tpR);
[map.aY, map.nY] = address(sz);              % addresses in y

sz = (segment(:,2)==tpS) + 3*(segment(:,2)==tpT) + ...
     (segment(:,2)==tpH) + 3*(segment(:,2)==tpR);
[map.aJ, map.nJ] = address(sz.*map.active);  % addresses in J

[a, n0] = address(sz);                       % data for fixed and movable
map.aP = a.*(segment(:,5)==1);               % p-addresses in x
map.aW = a.*(segment(:,5)==2);               % w-addresses in x

[a, n] = address(sz.*(segment(:,5)==1));     % data for movable only
map.nP = n;
map.nW = n0-n;

map.nX = map.nP + map.nW;                    % number of state variables


%---------- construct state info ----------------------------------------
tmp = [map.aP map.aW];                       % assemble x-addresses
z = zeros(map.nX,1);
f = false(map.nX,1);                         % initialize info structure
info = struct('type',z,'spatial',z);

for k = 2:nS                                 % loop over segments
   for tp = 1:2                              % loop over type (p,w)
      if tmp(k,tp)>0,                        % process if x-data exists
         adr = tmp(k,tp):(tmp(k,tp) + ...
            szj(segment(k,2)) - 1);          % make list of x-addresses

         info.type(adr) = tp;                % type (p=1,w=2)
         info.spatial(adr) = (segment(k,2)==tpH | ...
            segment(k,2)==tpR) + 1;          % spatial (trans=1,rot=2)
      end
   end
end

if nargout<=2,                               % matrices not needed
   return;
end


%------------ construct marices R, S, V ---------------------------------
datR = [datR datR/ratio].^2;                 % scale rotation component,
datS = [datS datS/ratio].^2;                 %  turn stdev into var
datV = [datV datV/ratio].^2;  

r = zeros(map.nX,1);                         % initialize diagonal of R
s = zeros(map.nX,1);                         % initialize diagonal of S
for tp = [2 1]                               % loop over w,p
   for sp = 1:2                              % loop over trans,rot
      ind = (info.type==tp & info.spatial==sp);
      r(ind) = dt * datR(3-tp+2*(sp-1));     % dynamics drift variance     
      s(ind) = datS(3-tp+2*(sp-1));          % initial state variance
   end
end
R = diag(r);                                 % make R from diagonal
S = diag(s);                                 % make S from diagonal

v = zeros(map.nR,1);                         % initialize diagonal of V
for k = 2:nS
   if segment(k,3)>0,                        % sensors only
      adr = map.aR(k):(map.aR(k) + ...
            szj(segment(k,3)) - 1);          % make list of r-addresses
      v(adr) = datV(1 + (segment(k,3)==tpH | ...
         segment(k,3)==tpR));                % residual variance
   end
end
V = diag(v);                                 % make V from diagonal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute addresses of segments' data in global array

function [a, n] = address(sz)

a = cumsum([1; sz]);
n = a(end)-1;
a = a(1:end-1) .* (sz>0);
