function [md,T]=detectstepsprep(varargin);
%
% [md,T]=detectstepsprep(mdata,markernames,Fs);
%
% Input
%    mdata        ->  Marker data. Matrix.
%    markernames  ->  Cell array of marker names.
%    Fs           ->  Sampling frequency.
     
% Returns the time series, and the approximate period to be used in the 
% automatic detection of cycles (steps). 

% Kjartan Halvorsen
% 2001-02-22
%
% Revisions
% 2002-12-10   Major revision. Changed to use popup menus instead
%              of input dialog

if (nargin==3) % Initial call. Construct a dialog
  mdata=varargin{1};
  mnames=varargin{2};
  Fs=varargin{3};

  h=detectstepsprep(mdata,mnames,Fs,[]);
  uiwait(h); % Resumed when user click 'ok' button

  mh=findobj(h,'Tag','markerPopup');
  markerind=get(mh,'Value');
  
  ch=findobj(h,'Tag','coordPopup');
  coordind=get(ch,'Value');

  dh=findobj(h,'Tag','durationEdit');
  if (~isempty(dh))
    duration=get(dh,'UserData'); % Should contain the step duration as a
				 % scalar
    T=Fs*duration;
  end 

  % Pick out the correct time series.
  md=mdata(:,(markerind-1)*3+coordind);

  % Close the dialog
  close(h);
    
elseif (nargin==4)  % Build the interface

  mdata=varargin{1};
  mnames=varargin{2};
  Fs=varargin{3};

  h0 = figure('Units','points', ...
	      'Color',[0.8 0.8 0.8], ...
	      'FileName','',...
	      'Name','Pick a time series', ...
	      'MenuBar','none', ...
	      'NumberTitle','off', ...
	      'PaperPosition',[10 10 140 140] , ...
	      'PaperUnits','points', ...
	      'Position',[200 200 165 175], ...
	      'Tag','Fig2', ...
	      'ToolBar','none');
  h1 = uicontrol('Parent',h0, ...
		 'Units','points', ...
		 'BackgroundColor',[0.8 0.8 0.8], ...
		 'ListboxTop',0, ...
		 'Position',[15 144 90 18], ...
		 'String','Use marker:', ...
		 'Style','text', ...
		 'HorizontalAlignment','left', ...
		 'Value',1);
  h1 = uicontrol('Parent',h0, ...
		 'Units','points', ...
		 'BackgroundColor',[0.5 0.6 0.69], ...
		 'ListboxTop',0, ...
		 'Min',1, ...
		 'Position',[15 126 135 17], ...
		 'String',mnames, ...
		 'Style','popupmenu', ...
		 'Tag','markerPopup', ...
		 'Value',1);
  h1 = uicontrol('Parent',h0, ...
		 'Units','points', ...
		 'BackgroundColor',[0.8 0.8 0.8], ...
		 'ListboxTop',0, ...
		 'Position',[15 101 90 18], ...
		 'String','Use coordinate:', ...
		 'Style','text', ...
		 'HorizontalAlignment','left', ...
		 'Value',1);
  h1 = uicontrol('Parent',h0, ...
		 'Units','points', ...
		 'BackgroundColor',[0.5 0.6 0.69], ...
		 'ListboxTop',0, ...
		 'Min',1, ...
		 'Position',[15 83 135 17], ...
		 'String',{'x','y','z'}, ...
		 'Style','popupmenu', ...
		 'Tag','coordPopup', ...
		 'Value',1);
  if (Fs>0)
    h1 = uicontrol('Parent',h0, ...
		   'Units','points', ...
		   'BackgroundColor',[0.8 0.8 0.8], ...
		   'ListboxTop',0, ...
		   'Position',[15 58 90 18], ...
		   'String','Duration of step (s):', ...
		   'Style','text', ...
		   'HorizontalAlignment','left', ...
		   'Value',1);
    h1 = uicontrol('Parent',h0, ...
		   'Units','points', ...
		   'BackgroundColor',[1 1 1], ...
		   'ListboxTop',0, ...
		   'Min',1, ...
		   'Position',[15 40 135 17], ...
		   'String','', ...
		   'Style','edit', ...
		   'Tag','durationEdit', ...
		   'HorizontalAlignment','left', ...
		   'Value',1);
  end
  h1 = uicontrol('Parent',h0, ...
		 'Units','points', ...
		 'BackgroundColor',[0.7 0.7 0.7], ...
		 'ListboxTop',0, ...
		 'Min',1, ...
		 'Position',[50 15 50 20], ...
		 'String','OK', ...
		 'callback','detectstepsprep(gcf)');
  md=h0; % Return a handle to the object
  
elseif (nargin==1) % Callback
  h=varargin{1};
  dh=findobj(h,'Tag','durationEdit');
  if (~isempty(dh))
    dstr=get(dh,'String'); % Should contain the step duration as a
			   % string
  
    duration=str2num(dstr);
    if ~isempty(duration)
      set(dh,'UserData',duration);
      uiresume
    end
  else
    uiresume
  end
  
elseif (nargin==0) 
  % Test this function

  % Create some data
  data=kron((1:9),ones(10,1));
  mnames={'m1','m2','m3'};
  Fs=240;
  
  [md,T]=detectstepsprep(data,mnames,Fs)
  
end
