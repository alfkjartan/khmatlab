function xv = qstate2vstate(xq, qi)
  %%  xv = qstate2vstate(xq, qi)
  %% Returns a state vector where the quaternion part is replaced by a
  %% 3x1 rotation vector. 
  %% Input 
  %%   xq   ->  state vector containing quaternion
  %%   qi   ->  index where quaternion starts
  %% Output
  %%   xv   <-  state vector witnout quaternion

  %% Kjartan Halvorsen 
  %% 2012-03-21

  if isempty(qi)
    xv = xq;
  else
    xv = zeros(size(xq,1)-1, size(xq,2));
    xv(1:qi-1,:) = xq(1:qi-1,:);
    xv(qi+3:end,:) = xq(qi+4:end,:);
    for i=1:size(xq,2)
      [w,th] = quaternion(xq(qi:qi+3, i));
      xv(qi:qi+2,i) = w'*th;
    end
  end


