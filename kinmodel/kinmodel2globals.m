function segm = kinmodel2globals(hm)
%  segm = kinmodel2globals(hm)
% Will set global variables, and return a set of segments, given the
% kinematic model hm (a struct).

% Kjartan halvorsen
% 2011-10-26

% Global variables. Ugly, but works if careful.
global GL_xi0 GL_g GL_gs GL_p0 GL_p GL_pnames GL_pdep

if (nargin == 0) % Unit test
  load unittestkinmodel % Loads hmtestmodel
  hmt = hmtestmodel;
  segmt = kinmodel2globals(hmt);

  %% Check
  thr = 1e-10;
  assert_equal("Correct dofs for root segment", length(segmt{1}.dofs), \
               length(hmt.twists{1}), thr);
  assert_equal("Correct twist for 1st dof", \
               GL_xi0(:,:,segmt{1}.dofs(1)), hmt.twists{1}{1}, thr);

  assert_equal("Correct dofs for second segment", length(segmt{2}.dofs), \
               length(hmt.twists{2}{1}), thr);
  %keyboard
  assert_equal("Correct twist for 1st dof of second segment", \
               GL_xi0(:,:,segmt{2}.dofs(1)), hmt.twists{2}{1}{1}, thr);
  ind = segmt{2}.markers(1);
  assert_equal("Correct initial 1st point  of second segment", \
               GL_p0(ind:ind+2), hmt.p0{2}{1}{1,2}, thr);
  markerinds = (segmt{2}.markers+2)/3;
  assert_equal("Correct dependencies for 2nd point  of second segment", \
               (1:8), GL_pdep{markerinds(2)} , thr);
  
else

  [segm, GL_xi0, GL_p0, GL_pnames, GL_pdep] = get_segments(hm.twists, hm.p0, \
                                                           0, 0, 0, []);

  GL_p0 = cat(1, GL_p0, ones(1,size(GL_p0, 2)));
  GL_g = zeros(size(GL_xi0));

  GL_p = GL_p0;

  GL_gs = zeros(4,4,length(hm.segm_names));
end


function [segments, xi0, p0, names, dependencies] = get_segments(tw, \
                                                                 p, \
                                                                 ind_s, \
                                                                 ind_xi, \
                                                                 ind_p, \
                                                                 predep)
  
  %% Recursively called 
  %%keyboard
    
  try
    
  if (size(p{1}, 2) > 1) 
    p0 = cat(2,p{1}{:,2});
  else
    p0 = [];
  end
  names = p{1}(:,1);
  xi0 = cat(3,tw{1}{:});
  dofs = size(xi0,3);
  segment.dofs = ind_xi + (1:dofs);
  segment.self = ind_s + 1;
  segment.markers = ind_p + (1:size(p0,2));
  segment.mnames = names;
  segment.parent = 0;
  segment.p0 = p0;

 catch
   keyboard
 end
  mydependencies = cat(2, predep, segment.dofs);
  dependencies = repmat({mydependencies}, length(segment.markers), 1);
  %%keyboard
  %markerind = (segment.markers + 2)/3;
  
  segments = cell(1,1);
  segments{1} = segment;

  ind_s = ind_s + 1;
  ind_xi = ind_xi + dofs;
  ind_p = ind_p + size(p0,2);

  if length(tw)>1 % Got branches
    depbrs = {};
    for br=2:length(p)
      [segmbr, xi0br, p0br, namesbr, depbr] = get_segments(tw{br},
                                                           p{br}, \
                                                           ind_s, ind_xi, \
                                                           ind_p, \
                                                           mydependencies); 
      try
      xi0 = cat(3, xi0, xi0br);
      p0 = cat(2, p0, p0br);
      names = cat(1, names, namesbr);
      depbrs = cat(1, depbrs, depbr);
      %%keyboard
      segmbr{1}.parent = segment.self;
      segments = cat(1, segments, segmbr);
      
      ind_s = ind_s + length(segmbr);
      ind_xi = ind_xi + size(xi0br,3);
      ind_p = ind_p + size(p0br,2);
      catch
        keyboard
      end
    end
    dependencies = cat(1,dependencies, depbrs); 
  end
endfunction






function dof = get_dofs(tw) 
  %% Recursively called 
  dof = length(tw{1});
  if length(tw)>1 % Got branches
    for br=1:length(tw{2})
      dof = dof + get_dofs(tw{2}(br));
    end
  endif
  endfunction

function xi0 = get_twists(tw) 
  %% Recursively called 
  xi0 = cat(3,tw{1}{:});
  
  for br=1:length(tw{2})
    xi0 = cat(3, xi0, get_twists(tw{2}{br}));
  end
endfunction

function [p0,names] = get_points(p)
  %% Recursively called 
  p0 = cat(1,p{1}{2,:});
  names = cat(1,p{1}{1,:});
  
  for br=1:length(p{2})
    [pbr,brn] = get_points(p{2}{br});
    p0 = cat(1, p0, pbr);
    names = cat(1, names, brn);
  end
endfunction

