function [Ms,Mstrue] = link_inertia(km)
%%  Ms = link_inertia(km)

if nargin == 0
   [Ms, Mstrue] = do_unit_test();
else
  if isstruct(km)
     %% Initial call
     II = km.inertia;
  else
      II = km;
  end
  try
  Ms = II{1};
  catch 
	keyboard
  end
  if length(II) > 1
    for i=2:length(II)
      Ms = cat(3, Ms, link_inertia(II{i}));
    end
  end
end


function [Ms, Mstrue] = do_unit_test()

km.inertia = {eye(6), {2*eye(6), {3*eye(6)}}};

Ms = link_inertia(km);

if (size(Ms,3) ~= 3)
   disp('Test 1 failed')
   disp('Expected'), 3
   disp('Found'), size(Ms,3)
else
   disp('Test 1 OK')
end

if ( Ms(6,6,3) ~= 3)
   disp('Test 2 failed')
   disp('Expected'), 3
   disp('Found'), Ms(:,:,3)
else
   disp('Test 2 OK')
end

km.inertia = {eye(6), ...
	      {2*eye(6), {3*eye(6)}},...
	      {4*eye(6), {5*eye(6), {6*eye(6)}}}};

Ms = link_inertia(km);

if (size(Ms,3) ~= 6)
   disp('Test 3 failed')
   disp('Expected'), 6
   disp('Found'), size(Ms,3)
else
   disp('Test 3 OK')
end

if ( Ms(6,6,3) ~= 3)
   disp('Test 4 failed')
   disp('Expected'), 3
   disp('Found'), Ms(:,:,3)
else
   disp('Test 4 OK')
end

if ( Ms(6,6,6) ~= 6)
   disp('Test 5 failed')
   disp('Expected'), 6
   disp('Found'), Ms(:,:,6)
else
   disp('Test 5 OK')
end
