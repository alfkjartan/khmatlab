function plotangles(states, names, a2plot, convertradians, unit)
%  plotangles(states,names) or
%  plotangles(states,names,angles2plot)
% 
% Plots the output from the EKF in separate windows.
%
% OBS: If the name of the angles does NOT contain a separated
% {x,y,z}, i.e. 'pelvis x' or 'pelvis_y', it is assumed that the
% value is in radians. The value is converted to degrees before
% plotting. 

% Input:
%     states      ->   matrices (nstates x nfr) with estimates
%                      state vectors.
%     names       ->   (optional) cell array (size nnstates) of strings.
%     a2plot      ->   (optional) nx2 cell
%                      array where the first column contains names
%                      of the angles to plot, and the second
%                      column contains the corresponding figure handles 
%     convertradians  ->  if 1, then the function will look for
%                      names that probably are angles, and convert
%                      from radians to degrees.
%     unit        ->   (optional) String appended to ylabel to
%                      indicate the unit of the values. Used mainly
%                      if used to plot joint angle contributions to
%                      end point velocity.
  

% Kjartan Halvorsen
% 2009-06-29   Based on plotmarkers from 2001-07-24


try
  
  [nsts,nfrs]=size(states);

  if (nargin < 5)
    unit = '';
  end
  
  if (nargin < 4)
    convertradians = 0;
  end
   
  if (nargin < 3)
  figs = zeros(nsts,1);
  for f=1:nsts
    figs(f) = figure('Name', names{f});
  end
  else
    if iscell(a2plot)
      states = cat(1, states,  zeros(1, nfrs));
      ind = [];
      nfound = 0;
      for i = 1:size(a2plot,1)
	[cmn,mind]=intersect(names,a2plot(i,1));
	if isempty(cmn)
	  ind=[ind nsts+1];
	else
	  nfound = nfound+1;
	  ind=[ind mind];
	end
      end
      states = states(ind,:);
      names = a2plot(:,1);
      figs = cat(1, a2plot{:,2});
      
    else
      figs = a2plot;
    end
  end

  T1=(1:nfrs);
  nsts = size(states,1);
  for i=1:min(length(figs),nsts)
      ystr = '';
    angle=names{i};
    figure(figs(i))
    if contains_str(angle, {'x','y','z', 'X', 'Y', 'Z'})
      conversion = 1;
      if ~isempty(unit)
	ystr = unit;
      end
    else
      if convertradians
	conversion = 180/pi;
	ystr = 'degrees';
      else
	conversion = 1;
	ystr= 'radians'
      end
      if ~isempty(unit)
	ystr = [unit];
      end
      
    end
    
    plot(T1,conversion*states(i,:));
   
    ylabel(ystr)
    xlabel('Frame');
    title(angle, 'Interpreter', 'none')
end

catch 
  keyboard
end

function has=contains_str(s, a) 

has=0;

for i=1:length(a)
  if ~isempty(strfind(s, [' ', a{i}]))
    has=1;
    return
  end
  if ~isempty(strfind(s, ['_', a{i}]))
    has=1;
    return
  end
end


