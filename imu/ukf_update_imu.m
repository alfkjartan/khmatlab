
function [M1,P1] = ukf_update_imu(M,P,Y,R, qi)

  %
  % Do transform and make the update
  %
  m = size(M,1); 

  %% create sigma points
  [X,Wm,Wc] = sigmaq(P,M,qi);

  %% propagate sigma points
  Z = observe_imu(X);


  [Xm, Zm, Px, Pz, Pxz] = sigma_cov_q(X, Z, qi, [], Wm, Wc);

  Pv = Pz + R;
  K = Pxz / Pv;
  v = Y-Zm;

  %% Update. Quaternion part is updated using quaternion multiplication
  M1 = qupdate(M, K, v, qi); 


  P1 = P - K * Pv * K'; % P - Pxz* Pv^-1 * Pv * Pv^-1 * Pzx 
		        % = P - Pxz * Pv^-1 * Pzx


  %keyboard