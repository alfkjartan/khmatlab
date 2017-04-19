function Aij = link_adjoint_trf(tws, th, i, j)
% Aij = link_adjoint_transformations(tws, th, i, j)
% Returns the adjoint transformations Ad (6 x 6 ) that transforms
% the twist of link j into the frame of  link i. 
%          Aij = \begin{cases} Ad^{-1}_(e^{\xi_{j+1}\theta_{j+1}}
%          \cdots e^{\xi_i \theta_i}), & i>j\\
%                   I_6, & i=j,\\
%                    0,  & i<j
% 
% Input
%   tw     ->  Array of twists \hat{\xi}_i\theta_i (4 x 4 x ndofs)
%              or (6 x ndofs) 
%   th     ->  Array of joint angles
%   i,j    ->  Indicies into the chain of twists.
% Output
%   Aij    <-  The adjoint transformations  (6 x 6) as defined above

%% Kjartan Halvorsen
% 2017-03-02

%% 
if nargin == 0
   do_unit_test();
else

  if i<j
      Aij = zeros(6,6);
      return 
  elseif i==j
      Aij = eye(6);
      return 
  else
      gprox = eye(4);
      for k = j+1:i
          if size(tws, 1) == 4
              gprox = gprox * expm(tws(:,:,k)*th(k));
          else
              gprox = gprox * expm(hat(tws(:,k))*th(k));
          end
      end
      [twslask, Aij] = adjoint_trf_inv(zeros(6,1), gprox);
      return
  end
end

end

%% 
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


Jsb0 = link_jacobians(tws, gsli, th);

twsflat = cat(3, hat([0;0;0;0;0;1]), ...
            hat([l1;0;0;0;0;1]), ...
            hat([l1+l2;0;0;0;0;1]), ...
            hat([0;0;1;0;0;0]));

gsliflat = cat(3, [I3 p1; 0 0 0 1], ...
                  [I3 p2; 0 0 0 1], ... 
                  [I3 p3; 0 0 0 1], ... 
                  [I3 p4; 0 0 0 1]);

for i=1:4
    Jsbi0 = Jsb0(:,:,i);
    Jsbi = zeros(6,4);
    for j=1:i
        Jsbi(:,j) = adjoint_trf_inv( link_adjoint_trf(twsflat, th, ...
                                                      i, j)*vee(twsflat(:,:,j)), ...
                                     gsliflat(:,:,i) );
    end
    tol = 1e-12;
    if norm(Jsbi0 - Jsbi) > tol
        disp(sprintf('Test %d failed', i))
        disp('Expected'), Jsbi0
        disp('Found'), Jsbi
        %keyboard
    else
        disp(sprintf('Test %d OK', i))
    end
end
end



