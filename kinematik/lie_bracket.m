function twn = lie_bracket(tw1, tw2)
% twn = lie_bracket(tw1, tw2)
% Lie bracket on se(3).  See section 3.3 Murray, Li, Sastry
%
% Arguments
%    tw1, tw2      --  twists, either 6x1 or 4x4

% Kjartan Halvorsen
% 2017-03-01

if nargin == 0
    do_unit_test()
else
    if size(tw1,1) == 6
        t1 = hat(tw1);
        t2 = hat(tw2);
        twn = vee(t1*t2 - t2*t1);
        return 
    else
        twn = tw1*tw2 - tw2*tw1;
        return 
    end
end

end

function do_unit_test()

    thr = 1e-12;
    
    tw1 = randn(6,1);
    tw2 = randn(6,1);
    tw3 = randn(6,1);
    
    % Test [tw1,tw2] = -[tw2,tw1]
    
    tw12= lie_bracket(tw1, tw2);
    tw21= lie_bracket(tw2, tw1);
    

    if norm(tw12+tw21) > thr
        disp('Test 1 failed.')
        disp('Expected  '), -tw12
        disp('Found  '), tw21
    else
        disp('Test 1 OK.')
    end
    
    % Test [tw1, [tw2, tw3]] + [tw2, [tw3, tw1]] + [tw3, [tw1, tw2]] = 0
    if norm( lie_bracket(tw1, lie_bracket(tw2, tw3)) ...
           + lie_bracket(tw2, lie_bracket(tw3, tw1)) ...
           + lie_bracket(tw3, lie_bracket(tw1, tw2)) ) > thr
        disp('Test 2 failed.')
        disp('Expected  0'), 
    else
        disp('Test 2 OK.')
    end

end
