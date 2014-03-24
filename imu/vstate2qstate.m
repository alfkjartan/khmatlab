function xq = vstate2qstate(xv, qi)
  %%  xq = vstate2qstate(xv, qi)
  %% Returns a state vector where the rotation vector part is replaced by a
  %% 4x1 quaternion
  %% Input 
  %%   xv   ->  state vector without quaternion
  %%   qi   ->  index where quaternion starts
  %% Output
  %%   xq   <-  state vector with quaternion

  %% Kjartan Halvorsen 
  %% 2012-03-21

if (nargin==0)
  do_unit_test();
else
  if isempty(qi)
    xq = xv;
  else
    w = xv(qi:qi+2);
    th = norm(w);
    w = w / th;
    q = quaternion(w, th);
    xq = cat(1, xv(1:qi-1), q', xv(qi+3:end));
  end
end

function do_unit_test()
  tol = 1e-12;

  disp('Unit test of vstate2qstate and qstate2vstate')

  qi = 9;
  x1 = randn(19,1);
  
  xq = vstate2qstate(x1,qi);
  xv = qstate2vstate(xq,qi);

  if (norm(xv-x1) > tol)
    disp('Test1: Failed')
    disp('Expected'), disp(x1)
    disp('Found'), disp(xv)
  else
    disp('Test1: OK')
  end


    

