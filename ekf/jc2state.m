function nkm=jc2state(km, initp, g_prox)
% function nkm=jc2state(km, initp, g_prox)
% Finds generalized coordinates that approximate the position of
% markers  state for the kinmodel that matches the given position of 
% the joint centers.
%
% Input
%    km       ->   kinmodel object
%    jc       ->   (tree) 3d coord of joint centers.
%    g_prox   ->   pose of proximal object.
% 
% Output
%    nkm      <-   the new kinmodel object.
%

% Kjartan Halvorsen
% 1999-09-20

jcc=jc{1};

g_prox;
invg=inv(g_prox);

brs=length(km.branches);

if (brs<2)
  Pjc=jcc;
  Pdjc=jc{2}{1};
  %The direction of the segment axis
  eax_c=invg*[(Pdjc-Pjc);0]; eax_c=eax_c(1:3)/norm(eax_c(1:3));
  eax_i=km.init_pose*[0;0;1;0]; eax_i=eax_i(1:3);
  % The rotations 
  th_rot=vects2var(km.jm,eax_i,eax_c);
  % the translation
  if (~isempty(km.j_sphere))
    Pjc_i=km.j_sphere(1:3);
  else
    Pjc_i=getCenter(km.jm);
  end
  Pjc=[Pjc;1]; Pjc=invg*Pjc;
  d=Pjc(1:3)-Pjc_i;
  [th_trl,flag]=trl2var(km.jm,d);
  km.state=[th_trl;th_rot];
else % Assume the first two branches give the lateral-medial
     % direction.
  Pjc=jcc;
  Pdjc1=jc{2}{1};
  Pdjc2=jc{3}{1};
  Pdjc=(Pdjc1+Pdjc2)/2;
  eax_c=Pdjc-Pjc; eax_c=eax_c(1:3)/norm(eax_c(1:3));
  elm_c=Pdjc1-Pdjc2; elm_c=elm_c(1:3)/norm(elm_c(1:3));
  eap_c=cross(elm_c,eax_c);

  g_total=[eap_c elm_c eax_c Pjc;0 0 0 1];
  % g_total=g_prox*km_current_pose*km_init_pose, which implies:
  g_c_p=invg*g_total*inv(km.init_pose);
  
  % Computing the state.
  [km.state,flag]=g2var(km.jm,g_c_p);
end

% compute the current pose of the kinmodel
km.current_pose=trf(km.jm,km.state);
g=g_prox*km.current_pose;

% Recursive call
for branch=1:brs
  km.branches{branch}=jc2state(km.branches{branch},jc{branch+1},g);
end

nkm=km;