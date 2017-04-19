function [eq, de_dw, deconj_dw] = qexp(w)
  %% [eq, de_dw, deconj_dw] = qexp(w)
  %% Computes a quaternian as the exponential of a vector w.

  %% Kjartan Halvorsen
  %% 2011-03-20

  if (nargin == 0)
    do_unit_test();
  else
    eq = zeros(4,1);
    
    wnorm = norm(w);
    if wnorm < 1e-10
      eq(4) = 1;
    else
      wn = w / wnorm;
      eq(1:3) = sin(wnorm/2)*wn;
      eq(4) = cos(wnorm/2);
    end

    if (nargout > 1)
      
      alpha = wnorm;
      %%alpha = 0;
      if (alpha < 1e-6)
	sa = 1/2 - alpha^2/48;
	ca = -1/24 + alpha^2/960;
	%%sa = 1;
	%%ca = -1/24;
      else
	sa = sin(alpha*0.5)/alpha;
	ca = 1/alpha^2 * ( 0.5*cos(alpha*0.5)-sa);
      end

      de_dw = cat(1, ca*w*w' + sa*eye(3), -sa*0.5*w');

      if (nargout > 2)
	Ic = -eye(4);
	Ic(4,4) = 1;
	deconj_dw = Ic*de_dw;
      end
    end % if nargout > 1
  end % if nargin == 0

function do_unit_test()

  w = randn(3,1);
  dw = 1e-6;
  I3 = eye(3);

  [ew,JJ,JC] = qexp(w);
  [ewc,JJC] = qexp(-w);

  
  de_dw = zeros(4,3);
  for i = 1:3
    de_dw(:,i) = (qexp(w+dw*I3(:,i)) - ew) / dw;
  end

  if ( norm(JJ - de_dw) > dw )
    disp('Test 1. Failed')
    cat(2, JJ, de_dw)
    keyboard
  else
    disp('Test 1. OK')
  end
 
  dec_dw = zeros(4,3);
  for i = 1:3
    dec_dw(:,i) = (qexp(-(w+dw*I3(:,i))) - ewc) / dw;
  end

  if ( norm(JC - dec_dw) > dw )
    disp('Test 2. Failed')
    cat(2, JC, dec_dw)
    keyboard
  else
    disp('Test 2. OK')
  end
 
  dt = 0.01;
  wdt = w*dt;
  [ew,JJ,JC] = qexp(wdt);
  [ewc,JJC] = qexp(-wdt);
  
  dewdt_dw = zeros(4,3);
  for i = 1:3
    dewdt_dw(:,i) = (qexp(dt*(w+dw*I3(:,i))) - ew) / dw;
  end

  dewdt_dw = zeros(4,3);
  for i = 1:3
    dewdt_dw(:,i) = (qexp(-dt*(w+dw*I3(:,i))) - ewc) / dw;
  end

  if ( norm(-dt*JJC - dewdt_dw) > dw )
    disp('Test 4. Failed')
    cat(2, -dt*JJC, dewdt_dw)
    keyboard
  else
    disp('Test 4. OK')
  end
 
  keyboard


  