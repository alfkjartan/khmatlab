function vv = qrot(v, q, invert)
%%  vv = qrot(v, q, invert)
%% Returns set of vectors or quaternions in v, rotated by given quaternion q.
%% 
%%   vv = q * v * q*
%
%% Kjartan Halvorsen
%% 2012-03-27

if (nargin == 0) 
  do_unit_test();
else

  if (nargin < 3)
    invert = 0;
  end

  nfr = size(v,2);

  vv = v;

  if (size(q,2) == nfr)
    if (size(v,1) == 3)
      for fr = 1:nfr
	if invert
	  qi = qinv(q(:,fr));
	else
	  qi = q(:,fr);
	end  
	vv(:,fr) = qtransv(v(:,fr), qi);
      endfor
    else
      for fr = 1:nfr
	if invert
	  qi = qinv(q(:,fr));
	else
	  qi = q(:,fr);
	end  
	vv(:,fr) = qtrans(v(:,fr), qi);
      endfor
    endif
  else
    if invert
      qi = qinv(q);
    else
      qi = q;
    end  
    if (size(v,1) == 3)
      for fr = 1:nfr
	vv(:,fr) = qtransv(v(:,fr), qi);
      endfor
    else
      for fr = 1:nfr
	vv(:,fr) = qtrans(v(:,fr), qi);
      endfor
    endif
  endif
endif


function do_unit_test()
  disp("Unit test of function qrot")
  
  qsb = quaternion([0;0;1], pi/2)'; % Rotates  vector in  body frame to spatial frame
  qbs = qinv(qsb); % Rotates  vector in spatial frame to body frame

  q1 = quaternion([1;0;0], pi/2)'; % Rotation of spatial frame

  q1b = qrot(q1,qbs); % rotation 1 in body frame. Should correspond to
		      % rotation about negative body a-axis
  q1bexpected = quaternion([0;-1;0], pi/2)';

  if (norm(q1b-q1bexpected) > 1e-12)
    disp('Test1: Failed')
    disp('Expected'), disp(q1bexpected)
    disp('Found'), disp(q1b)
  else
    disp('Test1: OK')
  end

  
  %keyboard
    
