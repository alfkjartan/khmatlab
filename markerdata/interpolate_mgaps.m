function nmd = interpolate_mgaps(md, order, window, maxgapsize, ...
				 mnames, fname)
%  nmd = interpolate_mgaps(md, order, window, maxgapsize, mnames)
% Function which interpolates multiple gaps in marker
% trajectories. The gaps are typically small but numerous in the
% specified window. A spline of specified order is applied. If a
% gap is larger than the specified max gapsize, it is left as it
% is.
% 
% Input
%    md          ->  marker data {attr, data}
%    order       ->  order of spline
%    window      ->  size of window
%    maxgapsize  ->  max size of gap
%    mnames      ->  optional cell array of marker names. If left
%                    out, all markers are processed.
% Output
%    nmd         <-  marker data

% kjartan Halvorsen
% 2007-04-12

if nargin == 0
  unittestme; % See below
  return
end

  
if nargin < 5
  mnames = getvalue(md{1}, 'MARKER_NAMES');
end

if nargin < 6
  fname = '';
else
  fname = fname(end-30:end);
end


nmd = md;

nfr = size(md{2}, 1);

halfwindow = floor(window/2);

for i=1:length(mnames)
  [traj,ind] = extractmarkers(md, mnames(i)); % nfr x 3 matrix
  ntraj = traj;
  
  % Find gaps. For each gap, center the window over the gap. Use
  % all data available within window to fit a polynom of desired
  % order. Present results
  
  plotflag = 0;
  
  if size(traj, 2) == 0
    continue
  end
  
  for comp=1:3
    nans = find(isnan(traj(:,comp)));
    index = 1;
    while (index < nfr & ~isempty(nans))
      % Find index of start and end of gap
      gapstart = nans(1);
      gapend = nans(1);
      gapindex = 1;
      if length(nans) > 1
	while (gapindex < length(nans))
	  gapindex = gapindex + 1;
	  if (gapend+1 == nans(gapindex))
	    gapend = gapend + 1;
	    if (gapindex == length(nans))
	      nans=[];
	    end
	  else
	    nans = nans(gapindex:end);
	    break
	  end
	end
      else
	nans = [];
      end
      
     
      gapcenter = floor((gapstart+gapend)/2);
      
      windowstart = max(1, gapcenter-halfwindow);
      windowend = min(nfr, gapcenter+halfwindow);
      
      win = (windowstart:windowend)';
      y = traj(win, comp);
      valueindex = find(~isnan(y));
      yy = y(valueindex);
      winwin = win(valueindex);
      p = polyfit(winwin,yy, order);
      
      ntraj(gapstart:gapend,comp) = polyval(p, gapstart:gapend);
      
      index = gapend+1;
    end
    
    if ~isempty(find(ntraj(:,comp) == 0))
      plotflag = 1;
    end
    
  end

  nmd{2}(:,ind) = ntraj;
  
  if plotflag
    fig=figure;
    clf
    subplot(311)
    plot((1:nfr), traj(:,1), 'b', 'LineWidth', 4);
    hold on
    plot((1:nfr), ntraj(:,1), 'r');
    subplot(312)
    plot((1:nfr), traj(:,2), 'b', 'LineWidth', 4);
    hold on
    plot((1:nfr), ntraj(:,2), 'r');
    subplot(313)
    plot((1:nfr), traj(:,3), 'b', 'LineWidth', 4);
    hold on
    plot((1:nfr), ntraj(:,3), 'r');
    subplot(311)
    title([fname, '  marker : ', mnames{i}])

    
    uiwait(fig)

  end
end

function unittestme
%keyboard

%pth = fileparts(mfilename('fullfile'));
%cd(pth)

%md = openmocapfile('',[pth '\..\..\..\projekt\gih\CoM_momentum\' ...
%		    'data\Gbrg0703\diskus070310\Veronica\' ...
%		    'kast0009_retrack.c3d']);

md = openmocapfile('', 'kast0009_retrack.c3d');
		    
md{2}(find(md{2}==0)) = NaN;

nmd = interpolate_mgaps(md, 3, 25, [], {'wrist_sup_r'});
		    