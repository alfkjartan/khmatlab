function [y,H]=observe_quick(x,yt,segments)
%  [y, H]=observe_quick(x, yt, segments)
% "Observation function" for a general mechanism
%  
% OBS: makes use of the following global variables:
%  GL_xi0 GL_g GL_gs GL_p0 GL_p GL_pnames
%
% Input
%    x        ->   the current state. (n x N) vector. Only the
%                  generalized coordinates, not time derivatives. 
%    yt       ->   the measurements.(3*m x 1) vector. Used for
%                  detecting missing markers.
%    segments ->   Cell array with segment structs
%                  dofs -> vector with indices into state vector
%                  markers -> vector with indices into marker vector
% Output
%    y     <-   the current marker positions. (3*m x 1) vector.
%    H     <-   Jacobian or the linearized observation matrix.
%               (3m x 2n) matrix. 

% Kjartan Halvorsen
% 2010-10-27
%
%

global GL_xi0 GL_g GL_gs GL_p0 GL_p GL_pnames GL_pdep


if (nargin>0)
  %% Compute the exponential of all the twists
  nn = size(x,1);
  for st=1:nn
    GL_g(:,:,st) = expr(GL_xi0(:,:,st)*x(st));
  end

  %% Compute the transformation for each segment and the position of the
  %% markers. 
  y = yt;

  for s=1:length(segments)
    segm = segments{s};

    par = segm.parent;
    if (par == 0)
      GL_gs(:,:,s) = eye(4);
    else
      GL_gs(:,:,s) = GL_gs(:,:,par);
    end

    for i=1:length(segm.dofs)
      GL_gs(:,:,s) = GL_gs(:,:,s)*GL_g(:,:,segm.dofs(i));
    end

    yp = GL_gs(:,:,s)*GL_p0(:,segm.markers);
 
    y = yp(1:3,:);
    y=y(:);

%    for m=1:length(segm.markers)
%      istart = segm.markers(m);
%      iend = istart+2;
%      y(istart:iend) = GL_gs(1:3,1:3,s)*GL_p0(istart:iend) + GL_gs(1:3,4,s);
%    end
  end

  return
  %% ------------------------------------------------
  %%  The Jacobian
  %% ------------------------------------------------


                                % transformations upto angle and after angle
  H = zeros(length(y), nn);

  stcomputed = zeros(nn,1);
  before=repmat(eye(4), [1 1 nn]);
                                %for k=1:nn-1
                                %  before(:,:,k+1)=before(:,:,k)*GL_g(:,:,k+1);
                                %  after(:,:,end-k)=GL_g(:,:,end-k+1)*after(:,:,end-k+1);
                                %end

  for m = 1:length(GL_pdep)
    after=repmat(eye(4), [1 1 nn]);
    mdeps = GL_pdep{m};
    indstart = (m-1)*3+1;
    indstop = m*3;
    p0 = GL_p0(indstart:indstop);
    for i=1:length(mdeps)
      st = mdeps(i);
      if ~stcomputed(st)
        for j=1:(i-1)
          before(:,:,st) = before(:,:,st)*GL_g(:,:,mdeps(j));
        end
        stcomputed(st)=1;
      end
      for j=i:length(mdeps)
        after(:,:,st) = after(:,:,st)*GL_g(:,:,mdeps(j));
      end
      dgdx = before(:,:,st)*GL_xi0(:,:,st)*after(:,:,st);
      H(indstart:indstop,st) = dgdx(1:3,1:3)*p0 + dgdx(1:3,4);
    end
  end
else % Unit test
  load unittestkinmodel % Loads hmtestmodel
  hmt = hmtestmodel;
  segmt = kinmodel2globals(hmt);

  x = zeros(8,1);
  x(1) = 1;

  [y,H] = observe_quick(x,zeros(3*7,1), segmt);
  npoints = length(y)/3;

  thr = 1e-10;
  dir1 = GL_xi0(1:3,4,1);
  %%keyboard
  assert_equal("First point correct", GL_p0(1:3)+dir1, y(1:3), thr);
  assert_equal("Jacobian for first state correct", repmat(dir1, \
                                                          npoints, 1), \
               H(:,1), thr);
  assert_equal("Last point correct", GL_p0(19:21)+dir1, y(19:21), thr);

  dt = 1e-6;
  x1 = zeros(8,1);
  x28 = x1;
  x28(8) = dt; 
  x24 = x1;
  x24(4) = dt;
  x27 = x1;
  x27(7) = dt;

  [y1,H] = observe_quick(x1,zeros(3*7,1), segmt);
  [y28,H] = observe_quick(x28,zeros(3*7,1), segmt);
  [y24,H] = observe_quick(x24,zeros(3*7,1), segmt);
  [y27,H] = observe_quick(x27,zeros(3*7,1), segmt);

  dydt8 = 1/dt*(y28-y1);
  dydt4 = 1/dt*(y24-y1);
  dydt7 = 1/dt*(y27-y1);
  assert_equal("Jacobian for fourth state correct", dydt4, H(:,4), 1e-4);
  assert_equal("Jacobian for seventh state correct", dydt7, H(:,7), 1e-4);
  assert_equal("Jacobian for last state correct", dydt8, H(:,8), 1e-4);
  keyboard
end
 