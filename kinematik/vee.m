function v=vee(vh)
%% For use with twists. See also hat.m
v = zeros(6,1);
v(1:3) = vh(1:3,4);
v(4) = vh(3,2);
v(5) = vh(1,3);
v(6) = vh(2,1);
