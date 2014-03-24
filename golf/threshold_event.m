function fr = threshold_event(traj, dir, thr, frames_after)
% Returns the first frame with threshold crossing

% Kjartan Halvorsen 
% 2010-03-24


dir = dir(:);

sig = traj*dir;

above_ind = find(sig > thr);


prev_ind = -1;
curr_ind = -1;
first_cross = nan;

n_above = -1;
for curr_ind = above_ind' 
  if n_above >= frames_after
    break;
  end
  
  if curr_ind == prev_ind+1
    n_above = n_above+1;
  else
    n_above = 1;
    first_cross = curr_ind;
  end
  
  prev_ind = curr_ind;
end

fr = first_cross;

  
