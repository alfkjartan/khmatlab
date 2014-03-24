function [msim, objd, mse] = simulate_and_plot (gm, states, dataframes, mdata, markers2plot, plotit)
	 
	 %% usage: [msim, objsim, mse] = simulate_and_plot (gm, states, dataframes, mdata, markers2plot, plotit)
	 %% 
	 %% 
	 
%% Kjartan Halvorsen
%% 2013-10-09

debug = 1; % Will compute end point trajectory with forward_map.m and return as third
           % output argument

[nsts, nfrs] = size(states);
[msim, simnames, objd, objnames] = sim_model(gm, states(1:nsts/2,:), 'objectcenter');

if debug
   mse = zeros(size(objd));
   nbranches = size(objd,2)/3;
   g0 = gm.object_frame;
   for i=1:nfrs
     gsti = forward_map(gm.twists, g0, states(1:nsts/2, i));
     for br = 1:nbranches
	 mse(i, (br-1)*3+1:br*3) = gsti(1:3, 4, br)';
     end
   end
end


if plotit
  % Plot results
  y_observations=extractmarkers(mdata, simnames);
  y_observations(find(y_observations==0)) = NaN;
  msee = plotmarkers(y_observations(find(dataframes),:), ...
 		    simnames,...
 		    msim, simnames, markers2plot, {'data', 'model'});
  
  if ~debug
     mse = msee;
  end
end




