
function [M1,P1] = ukf_update_d(M,P,Y,R)

  %
  % Do transform and make the update
  %
  m = size(M,1); 

  %% create sigma points
  [X,Wm,Wc] = sigma(P,M,0);

  %% propagate sigma points
  %Z = observe_imu(X);
  Z = X(7:9,:);


  [Xm, Zm, Px, Pz, Pxz] = sigmacov(X, Z,  0, 0, Wm, Wc);

  Pv = Pz + R;
  K = Pxz / Pv;
  v = Y-Zm;

  %% Update. Quaternion part is updated using quaternion multiplication
  M1 = M + K*v;


  P1 = P - K * Pv * K'; % P - Pxz* Pv^-1 * Pv * Pv^-1 * Pzx 
		        % = P - Pxz * Pv^-1 * Pzx


  %keyboard