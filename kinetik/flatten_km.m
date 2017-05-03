function [tws, g0, Mb, gcnames] = flatten_km(km)
% 
% Returns a flattened version of the kinematic model km. By flatten
% means that the kinematic chain is represented as a series of
% links with 1-dof joints inbetween. Some of these links have zero
% length and zero mass.
%
% Input
%    km   ->  kinematic model with fields
%             twists -> nested array of twists
%             g0     -> nested array of transformations from local
%                       to spatial frame 
%             mass   -> nested array of masses
%             inertia -> nested array of inertia matrices
%             CoM -> nested array of CoM position
% Output
%    tws  <-  Array of twists (6 x ndofs)
%    g0   <-  Array of g0s (4 x 4 x ndofs)
%    Mb   <-  Array of local inertia matrices (6 x 6 x ndofs)
%    gcnames <- list of names of the dofs

%% Kjartan Halvorsen
%  2017-03-20
    
if nargin == 0
    do_unit_test()
else
    linkinds = link_indices(km.twists); 
    tws = vee(unravel(km.twists));
    ndofs = size(tws, 2);

    g0 = repmat(eye(4), [1,1,ndofs]);
    g0_tmp = unravel(km.g0);
    g0(:,:,linkinds) = g0_tmp;

    Mb = zeros(6,6,ndofs);
    Mb_tmp = unravel(km.inertia);
    Mb(:,:,linkinds) = Mb_tmp;

    gcnames = km.gcnames;
    
end

end

function i = link_indices(vnested, proxdofs)
% Will go through the nested array vnested and determine indices
% into the flattened array where non-zero mass links are
    if nargin == 1
        proxdofs = 0;
    end
    mydofs = length(vnested{1});

    i = proxdofs + mydofs;
    
    if length(vnested) > 1
        i = cat(1, i, link_indices(vnested{2}, proxdofs+mydofs));
    end
end

function v = unravel(vnested)
% Will unravel the nested arrays in vnested
    vv = vnested{1};
    if iscell(vv)
        ndim = size(vv{1}, 2);
        
        if ndim > 1
            % Matrix so stack along 3rd dimension
            v = cat(3, vv{:});
        else
            % Vector so stack along 2nd dimension
            v = cat(2, vv{:});
        end
    else
        v = vv;
        ndim = size(v,2);
    end
    if length(vnested) > 1
        if ndim > 1
            v = cat(3, v, unravel(vnested{2}));
        else
            v = cat(2, v, unravel(vnested{2}));
        end
    end
end


function do_unit_test()
    
    km = scara_robot_model();
    
    [tws, g0, Mb] = flatten_km(km);
    
    if size(tws, 1) ~= 6 
        disp('Test 1 failed')
        disp('Expected twists to be length 6')
        disp('Found'), size(tws, 1)
        tws
    else
        disp('Test 1 OK')
    end

    if size(tws, 2) ~= 4 
        disp('Test 2 failed')
        disp('Expected ndofs to be 4')
        disp('Found'), size(tws, 2)
        tws
    else
        disp('Test 2 OK')
    end

    if size(g0, 1) ~= 4 
        disp('Test 3 failed')
        disp('Expected transformation')
        disp('Found'), g0
        tws
    else
        disp('Test 3 OK')
    end
    if size(g0, 2) ~= 4 
        disp('Test 4 failed')
        disp('Expected transformation')
        disp('Found'), g0
        tws
    else
        disp('Test 4 OK')
    end
    if size(g0, 3) ~= 4 
        disp('Test 5 failed')
        disp('Expected transformation')
        disp('Found'), g0
        tws
    else
        disp('Test 5 OK')
    end

end