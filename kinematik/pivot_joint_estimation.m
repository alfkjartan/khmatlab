function [jc, U, S, V] = pivot_joint_estimation(refdata, mdata, ...
						proxmarkers, ...
						distmarkers) 
%  [jc, U, S, V] = pivot_joint_estimation(refdata, mdata, ...
%						proxmarkers, ...
%						distmarkers) 
% 
% Implementation of the symmetric pivot method for estimating the
% center of rotation.
%
% Input
%    refdata   ->  {attr, data} Mocap data with reference position
%                  of the markers
%    mdata     ->  {attr, data} Mocap data with movement
%    proxmarkers ->  cell array with marker names
%    distmarkers ->  cell array with marker names
% Output
%    jc        <-  the estimated joint center. jc is [jcd;jcp],
%                  with jcp and jcd denoting the position of the joint
%                  center in the proximal and distal, respectively,
%                  segment's reference system.
%    U,S,V     <-  the svd of the left hand side matrix

% Kjartan Halvorsen
% 2006-02-15

usedecimate=0;

nfr = size(mdata{2},1);
ncol = size(mdata{2},2);

if (nfr>1000) % Decimate data (downsample) to avoid memory problem
  samplefreq = str2double(getvalue(mdata{1}, 'FREQUENCY'));
  resfreq = 2; %Hz
  fff =  fix(samplefreq/resfreq);

  if usedecimate
    factorvec = factor(fff);
    for ff=factorvec
      mdc = decimate(mdata{2}(:,1), ff);
      mdcopy = zeros(length(mdc), ncol);
      mdcopy(:,1) = mdc;
      for c=2:ncol
	mdcopy(:,c) = decimate(mdata{2}(:,c), ff);
      end
      mdata{2} = mdcopy;
    end
  else
    mdata{2} = downsample(mdata{2}, fff);
  end

end % if (nfr>1000)

%Tprox = getMotion(refdata, mdata, proxmarkers);
%Tdist = getMotion(refdata, mdata, distmarkers);
Tprox = getMotion( mdata, proxmarkers);
Tdist = getMotion( mdata, distmarkers);

T0prox = getMotion( refdata, mdata, proxmarkers);
Tp0 = reshape(T0prox(1,:), 4,4);
Tp0inv = ginv(Tp0);

nfr = size(Tprox,1);

A = zeros(nfr*3, 6);
b = zeros(nfr*3,1);

for i=1:nfr
  Tpi = reshape(Tprox(i,:), 4,4);
  Tdi = reshape(Tdist(i,:), 4,4);

  if (~Tdi(1:3, 1:3) | ~Tpi(1:3, 1:3))
    % markers are missing. Skip this frame
    else
      A((i-1)*3+1:(i*3), :) = cat(2, Tdi(1:3,1:3), -Tpi(1:3, 1:3));
      b((i-1)*3+1:(i*3)) = Tpi(1:3,4)-Tdi(1:3,4);
  end
end
[U,S,V] = svd(A);

Ainv = generalized_inverse(U,S,V);
jc = Ainv*b;

jc1 = Tp0inv*[jc(1:3);1];
jc2 = Tp0inv*[jc(4:6);1];

jc = cat(1, jc1(1:3), jc2(1:3));

V1 = Tp0inv*[V(1:3,6);0];
V2 = Tp0inv*[V(4:6,6);0];

V(:,6) = cat(1, V1(1:3), V2(1:3));


%keyboard