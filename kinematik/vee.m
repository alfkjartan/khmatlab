function v = vee(vh)
% For use with twists. See also hat.m
ntws = size(vh, 3);

v = zeros(6,ntws);

for i=1:ntws
    v(1:3,i) = vh(1:3,4,i);
    v(4,i) = vh(3,2,i);
    v(5,i) = vh(1,3,i);
    v(6,i) = vh(2,1,i);
end
