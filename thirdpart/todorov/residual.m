%[r, predict, J] = residual(y, x, segment, map)
%
% Compute residual and Jacobian of kinematic model
%
% y: measurement column vector, NaN indicates missing data
% x: state column vector
%
% segment: nS x 5 array of segment descriptors (rows) with elements:
%   1. index of parent segment (1- parent is root)
%   2. joint type: 1- S(lide); 2- T(ranslation); 3- H(inge), 4- R(otation)
%   3. sensor type, same convention as joint type; 0- no sensor
%   4. axis of S or H joint: 1- X_axis, 2- Y_axis, 3-Z_axis
%   5. joint category:  0- none, 1- movable, 2- fixed
%
% map: structure (usually created by prepare.m) with fields:
%   active: nS x 1 logical array specifying joints to be included in J
%   aR, aY, aJ: nS x 1 integer arrays specifying addresses of the
%      joints' data in the vectors r,y and matrix J (column index)
%   aP, aW: nS x 1 integer arrays specifying addresses of the
%      joints' p,w data in the state vector x; 0- joint has no data
%   nR, nY, nJ, nP, nW: total amout of data in each category
%   nX: size of state vector; nX = nP+nW
%
% r: residual vector, NaN when corresponding data are missing in y
% predict: predicted measurement vector, same format as y
% J: Jacobian matrix (d r / d active_joints) 

% Copyright (C) Emanuel Todorov, 2006-2007

function [r, predict, J] = residual(y, x, segment, map)

tpS = 1;                                     % type Slide
tpT = 2;                                     % type Translation
tpH = 3;                                     % type Hinge
tpR = 4;                                     % type Rotation

nS = size(segment,1);                        % number of segments

axis = zeros(3,nS);                          % allocate S/H joint axis
axis(1,segment(:,4)==1) = 1;                 % set to X/Y/Z unit vector
axis(2,segment(:,4)==2) = 1;
axis(3,segment(:,4)==3) = 1;

aX = map.aP + map.aW;                        % p,w-addresses in x


%------------ compute forward kinematics, residuals and predictions ------
r = NaN*zeros(map.nR,1);                     % allocate residuals
predict = NaN*zeros(map.nR,1);               % allocate predictions

frmPos = zeros(3,nS);                        % allocate frame positions
frmQuat = zeros(4,nS);                       % allocate frame quaternions
frmRot = zeros(3,3,nS);                      % allocate frame rotations

frmRot(:,:,1) = eye(3);                      % initialize root frame
frmQuat(:,1) = [1;0;0;0];

for k=2:nS                                   % loop over non-root segments
   par = segment(k,1);                       % get segment info
   jnt = segment(k,2);
   stype = segment(k,3);
   dat = map.aP(k)+map.aW(k);
   res = map.aR(k);
   sen = map.aY(k);

   frmPos(:,k) = frmPos(:,par);              % copy frame from parent
   frmQuat(:,k) = frmQuat(:,par);
   frmRot(:,:,k) = frmRot(:,:,par);

   %-------- compute forward kinematics
   switch jnt,
      case tpS,                              % Slide joint
         frmPos(:,k) = frmPos(:,par) + frmRot(:,:,par)*(x(dat)*axis(:,k));

      case tpT,                              % Translation joint
         frmPos(:,k) = frmPos(:,par) + frmRot(:,:,par)*x(dat:dat+2);

      case tpH,                              % Hinge joint
         frmQuat(:,k)=QuatMult(frmQuat(:,par),Axis2Quat(x(dat)*axis(:,k)));
         frmRot(:,:,k) = Quat2Matrix(frmQuat(:,k));

      case tpR,                              % Rotation joint
         frmQuat(:,k) = QuatMult(frmQuat(:,par), Axis2Quat(x(dat:dat+2)));
         frmRot(:,:,k) = Quat2Matrix(frmQuat(:,k));
         
      otherwise
         error(sprintf('invalid joint type: %d', jnt));
   end

   %-------- compute residuals and predictions
   if stype && ~isnan(y(sen)),               % only available sensors
      switch stype,
         case tpT,                           % Translation sensor
            r(res:res+2) = y(sen:sen+2) - frmPos(:,k);
            predict(sen:sen+2) = frmPos(:,k);

         case tpH,                           % Hinge sensor
            r(res) = y(sen) - x(dat);
            predict(sen) = x(dat);

         case tpR,                           % Rotation sensor
            qdif = QuatMult(QuatNeg(frmQuat(:,k)), y(sen:sen+3));
            r(res:res+2) = qdif(2:4);
            predict(sen:sen+3) = frmQuat(:,k);
            
         otherwise
            error(sprintf('invalid sensor type: %d', jnt));
      end
   end
end


%--------------- compute Jacobian ----------------------------------------
if nargout<3,
   return;                                   % exit if J not needed
end

J = zeros(map.nR,map.nJ);                    % allocate Jacobian
J(isnan(r),:) = NaN;                         % NaN for missing data

for s = 2:nS                                 % loop over non-root segments  
   stype = segment(s,3);                     % get sensor info
   res = map.aR(s);
   sen = map.aY(s);

   if ~stype || isnan(y(sen)),               % skip unavailable sensors
      continue;                              
   end

   if stype==tpH,                            % fast processing of H sensor
      if segment(s,2)==tpH,
         if map.active(s),
            J(res,map.aJ(s)) = -1;           % assume offset=0, gain=1
         end
         continue;
      else
         error('goniometer can only be attached to hinge joint');
      end
   end

   k = s;
   while k>1,                                % loop over ancestor joints
      if ~map.active(k),
         k = segment(k,1);                   % move to parent
         continue;                           % skip inactive joint
      end

      par = segment(k,1);                    % get joint info
      jnt = segment(k,2);
      dat = map.aP(k)+map.aW(k);
      jac = map.aJ(k);

      % ----------- compute derivatives of joint quaternion
      if jnt==tpH,                           % Hinge joint
         dq = [-sin(x(dat)/2); cos(x(dat)/2)*axis(:,k)]/2;
         dq1 = [dq(1); -dq(2:4)];            % conjugate

      elseif jnt==tpR,                       % Rotation joint
         dq = QuatDeriv(x(dat:dat+2));
         dq1 = [dq(1,:); -dq(2:4,:)];        % conjugate
      end

      % ----------- compute J for all sensor/joint combinations
      if stype==tpT,                         % Translation sensor
         switch jnt,
            case tpS,                        % Slide joint
               J(res:res+2,jac) = -frmRot(:,:,par)*axis(:,k);

            case tpT,                        % Translation joint
               J(res:res+2,jac:jac+2) = -frmRot(:,:,par);

            case {tpH,tpR},                  % Hinge/Rotation joint
               qV = [0; frmRot(:,:,k)'*(frmPos(:,s)-frmPos(:,k))];
               if jnt==tpH,
                  qA = Axis2Quat(x(dat)*axis(:,k));
               else
                  qA = Axis2Quat(x(dat:dat+2));
               end

               for t=1:size(dq,2)            % process all columns of dq
                  tmp = QuatMult(frmQuat(:,par), QuatMult(...
                         (QuatMult(dq(:,t),QuatMult(qV,QuatNeg(qA))) + ...
                          QuatMult(qA,QuatMult(qV,dq1(:,t)))), ...
                         QuatNeg(frmQuat(:,par))));

                  J(res:res+2,jac+t-1) = -tmp(2:4);
               end               
         end

      elseif stype==tpR,                     % Rotation sensor
         if jnt==tpH || jnt==tpR,            % Hinge/Rotation joint
            qT = QuatMult(QuatNeg(frmQuat(:,s)), frmQuat(:,k));
            qP = QuatMult(QuatNeg(frmQuat(:,par)), y(sen:sen+3));

            for t=1:size(dq1,2)              % process all columns of dq1
               tmp = QuatMult(qT, QuatMult(dq1(:,t),qP));
               J(res:res+2,jac+t-1) = tmp(2:4);
            end
         end
      end

      k = par;                               % move to parent
   end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiply two quaternions

function q = QuatMult( q1, q2 )

q = [q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);...
     q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3);...
     q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);...
     q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)];

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert unit quaternion to 3x3 orthonormal matrix

function mat = Quat2Matrix( q )

mat = [q(1)^2+q(2)^2-q(3)^2-q(4)^2 ...
       2*(q(2)*q(3)-q(1)*q(4)) ...
       2*(q(2)*q(4)+q(1)*q(3)); ...
       
       2*(q(2)*q(3)+q(1)*q(4)) ...
       q(1)^2-q(2)^2+q(3)^2-q(4)^2 ...
       2*(q(3)*q(4)-q(1)*q(2)); ...
          
       2*(q(2)*q(4)-q(1)*q(3)) ...
       2*(q(3)*q(4)+q(1)*q(2)) ...
       q(1)^2-q(2)^2-q(3)^2+q(4)^2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form conjugate quaternion

function qn = QuatNeg( q )

qn = [q(1); -q(2); -q(3); -q(4)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form unit quaterion from 3D vector whose length is the rotation angle

function q = Axis2Quat( axis )

angle = sqrt(axis'*axis);
if angle>0,
    q = [cos(angle/2); sin(angle/2)*axis/angle];
else
    q = [1; 0; 0; 0];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Jacobian matrix of exponential map

function J = QuatDeriv( axis )

ang = norm(axis);

if ang<1E-5,
    sa = 1/2 - ang^2/48;
    ca = -1/24 + ang^2/960;
else
    sa = sin(ang/2)/ang;
    ca = (cos(ang/2)/2 - sa)/ang^2;
end

J = [-sa/2; ca*axis] * axis';
J(2:4,:) = J(2:4,:) + eye(3)*sa;
