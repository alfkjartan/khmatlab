function [gg, pb] = paddel_motion(gs, md, pks, bf_fcn, bf_markers)
%%  [gg, pb] = paddel_motion(gs, md, pks, bf_fcn, bf_markers)
%% Computes the movement of the paddel from motion capture data. 
%%
%% Input
%%   gs                ->  transformation from first frame of data to
%%                         all subsequent frames
%%   md                ->  marker data
%%   pks               ->  list of indices where cycle starts
%%   bf_fcn            ->  name of function which computes local
%%                         coordinate frame
%%   bf_markers        ->  struct used by bf_fcn

%% Kjartan Halvorsen
%% 2012-05-10

nfrs = size(gs,1);
gg = zeros(4,4,nfrs);
pb = zeros(nfrs, 4);

for i=1:length(pks)
  p = pks(i);

  %% Compute the transformation between local frame at pos p and frame 1
  [gsb, blade] = feval(bf_fcn, md, bf_markers, p); 
  gbs = ginv(gsb);

  %blade
  
  gp0 = reshape(gs(p,:), 4, 4);
  g0p = ginv(gp0);


  %% Transform all transformations to the static body frame between peak
  %% p and next peak
  if (i < length(pks))
    stopfr = pks(i+1)-1;
  else
    stopfr = size(gs,1);
  end

  %keyboard

  for fr = p:stopfr
    gfr0 = reshape(gs(fr,:), 4, 4);
    gg(:,:,fr) = gbs*gfr0*g0p*gsb;
    pb(fr,:) = (gg(:,:,fr)*blade)';
  end
end

pb(:,4) = [];


  
