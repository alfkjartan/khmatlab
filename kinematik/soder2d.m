function [T,res]=soder2d(data)
% function [T,res]=soder2d(data)
%
% Computes 3 x 3 transformation matrix for 2d data. Modification of
% Söderkvist & Wedin's method for 3d.

% Kjartan Halvorsen
% 2004-04-29

if (nargin == 1)
  m2 = size(data,2);

  if mod(m2,2)
    disp('ERROR: input has to be multiple of 2 (XY coordinates)'); return
  end

  A=[reshape(data(1,:)',2,m2/2)]';  
  B=[reshape(data(2,:)',2,m2/2)]';  


  % Checking for NaNs and also checking if still 3 pts left and if not
  % T=[NaN...];
  cut=[0];
  qA=isnan(A); qB=isnan(B); qAB=[qA,qB];
  qsum=sum(qAB'); cut=find(qsum~=0);
  A([cut],:)=[]; B([cut],:)=[];
  if size(A,1)<2,
    T = repmat(NaN,3,3);
  end


  Amean=mean(A)'; Bmean=mean(B)';


  for i=1:size(A,1)-size(cut,2),
    Ai(:,i)=[A(i,:)-Amean']';
    Bi(:,i)=[B(i,:)-Bmean']';
  end


  C=Bi*Ai';
  [P,T,Q]=svd(C);

  R=P*diag([1 det(P*Q')])*Q';

  d=Bmean-R*(Amean);

  T=[R d; 0 0 1];

  % Calculating the norm of residuals
  A = cat(1, A', ones(1,size(A,1)));
  B=B';
  Bcalc=T*A; 
  Diff=B-Bcalc(1:2,:); Diffsquare=Diff.^2;
  %DOF=2*(number of points)-3 unknowns (dx,dy,alpha):
  DOF=max(size(B,1)*size(B,2)-3,1);
  res=[sum(Diffsquare(:))/DOF].^0.5;

else % Unit test
  disp('Unit test of soder2d function')
  
  % Tolerance
  tol = repmat(1e-12,3,3);
  
  % Set up test case.
  tlngth = 400;

  % Translation
  d(1:3,1,:) = cat(1, 2*rand(2, tlngth) - 1, ones(1,tlngth));
  
  % Rotation
  phi = linspace(0,5*pi/6,tlngth);
  cp(1,1,:) = cos(phi);
  sp(1,1,:) = sin(phi);
  zs = zeros(1,2,tlngth);
  R = [ [cp -sp; sp cp; zs] d ]; 

  % Move some points
  p0 = [1 -1 0 0
	0 0 1 -1
	1 1 1 1];
  
  p0rad = p0(1:2,:);
  p0rad = (p0rad(:))';
  
  ok=1;
  
  for i=1:tlngth
    p = R(:,:,i)*p0;
    
    T = soder2d(cat(1,p0rad, reshape(p(1:2,:),1,8)));
    
    oki = min( min( abs( T - R(:,:,i) ) < tol));
    
    if ~oki
      disp(['Assertion failed'])
      expected = R(:,:,i)
      found = T
    end

    ok = ok & oki;
  end

  T = ok;
end




