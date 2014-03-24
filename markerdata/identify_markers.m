function alilist=identify_markers(varargin)
%  alilist=identify_markers(markerid,im,markernames)
% Lets the user define the identity of each marker in the set.
%

% Some code borrowed from Gait/cb_gait_markers_setup.m
  
% Kjartan Halvorsen
% 2002-12-16

if (nargin==0) % Testing this function
  testfunction
  return
end

if (nargin>2) % Initial call, set up the dialog

  % check the arguments
  if (checkargs(varargin)) % Returns 1 if not ok
    return
  end
  
  markerid=varargin{1};
  nmid=size(markerid,1);
  im=varargin{2};
  mnames=varargin{3};
  
  startindices=ones(nmid,1);
  
  if nargin==4 % Defaults for some markers are given
    defaults=varargin{4};
    for i=1:nmid
      [inms,ind]=intersect(mnames,defaults{i,2});
      if ~isempty(ind)
	startindices(i)=ind;
      end
    end
  end
  
  
  % Check the markerid array for the most closely spaced list
  % positions. Set image size accordingly.
  minleft=realmax;
  minright=realmax;
  for i=1:nmid
    for j=i+1:nmid
      if (markerid{i,2}(1)==markerid{j,2}(1))
	if (markerid{i,2}(1)==0) % Left edge
	  if (abs(markerid{i,2}(2)-markerid{j,2}(2))<minleft)
	    minleft=abs(markerid{i,2}(2)-markerid{j,2}(2));
	  end
	else
	  if (abs(markerid{i,2}(2)-markerid{j,2}(2))<minright)
	    minright=abs(markerid{i,2}(2)-markerid{j,2}(2));
	  end
	end
      end
    end
  end
  
  minspace=min(minright,minleft);
  if (minspace<1e-8)
    errordlg({'Error in function identify_markers',...
	      'Markerpositions on image too close.'});
    return
  end

  [h,w,rgb]=size(im);
  image_height=fix(27/minspace);
  image_width=image_height*w/h; %get right image aspect ratio

  pu_width=100; %popup width
  pu_height=15; %popup height
  align_left=20;
  align_right=align_left*2+pu_width+image_width;
  align_diff=align_left+pu_width+image_width;

  frame1_height=image_height+20;

  fig_height=frame1_height+70;
  fig_width=align_right+pu_width+20;

  frame1_y=fig_height-frame1_height-10;

  %draw the figure
  markers_fig=figure('visible','off');
  set(markers_fig,'tag','markers_fig')
  %  set(markers_fig,'position',get_screen_pos(fig_width,fig_height)) 
  set(markers_fig,'position',[100 100 fig_width fig_height]) 
  set(markers_fig,'color',[0.8 0.8 0.8])
  set(markers_fig,'menubar','none')
  set(markers_fig,'NumberTitle','off');
  set(markers_fig,'Name','Marker Setup');
  %  set(markers_fig,'windowstyle','modal')

  row_pos=frame1_y+265;

  %draw the marker man
  ax=axes;
  set(ax,'tag','ax')
  set(ax,'unit','pixels')
  ax_w=fig_width-10;
  set(ax,'position',...
	 [5 fig_height-frame1_height-10 ax_w frame1_height])
  imHand=image(im);
  xdata=get(imHand,'XData');
  ydata=get(imHand,'YData');
  set(ax,'xlim',[0 ax_w])
  set(ax,'ylim',[0 frame1_height])
  set(imHand,'XData',[(ax_w-image_width)/2 (ax_w+image_width)/2])
  set(imHand,'YData',[10 image_height+10])

  %input frame
  l1=line([0 0],[0 frame1_height-0.9]);
  set(l1,'color',[0 0 0])
  l2=line([0 ax_w],[0 0]);
  set(l2,'color',[0 0 0])
  l3=line([ax_w ax_w],[0 frame1_height]);
  set(l3,'color',[0 0 0])
  l4=line([0 ax_w],[frame1_height-0.99 frame1_height-0.99]);
  set(l4,'color',[0 0 0])
 
  set(ax,'visible','off')

  %draw popup menus for the markers
  for i=1:nmid
    halign=align_left+markerid{i,2}(1)*align_diff;
    vpos=60+fix(markerid{i,2}(2)*image_height);
    hh=uicontrol('style','popupmenu',...
		 'tag',markerid{i,1},...
		 'string',mnames,...
		 'value',startindices(i),...
		 'position',[halign vpos pu_width pu_height],...
		 'enable','on');
  end
  

  cancel_btn=uicontrol('style','pushbutton');
  set(cancel_btn,'tag','cancel_btn')
  set(cancel_btn,'position',[fig_width-70-70-75 15 60 25])
  set(cancel_btn,'string','Cancel')
  set(cancel_btn,'callback','identify_markers(gcf,''cancel'')')

  back_btn=uicontrol('style','pushbutton');
  set(back_btn,'tag','back_btn')
  set(back_btn,'position',[fig_width-70-75 15 60 25])
  set(back_btn,'string','Ok')
  set(back_btn,'callback','identify_markers(gcf,''ok'')')

  stath=uicontrol('style','text');
  set(stath,'tag','status')
  set(stath,'BackgroundColor',[0.8 0.8 0.8])
  set(stath,'position',[15 15 fig_width-70-70-75-25 25])

  set(markers_fig,'visible','on')

  set(markers_fig,'UserData',markerid(:,1));
  
  drawnow

  uiwait(markers_fig);
  
  % Returns here when finished.
  ud=get(markers_fig,'UserData');
  if (~isempty(ud))
    if (ischar(ud))
      if (strcmp(ud,'cancelled'))
	alilist={};
      end
    else
      % Return list
      alilist=ud;
    end
  else
    alilist={};
  end
  
  close(markers_fig)
  
elseif (nargin==2) % Callback
  figh=varargin{1};
  action=varargin{2};

  stath=findobj(figh,'tag','status');
  set(stath,'String','');

  switch action
   case 'cancel'
    set(figh,'UserData','cancelled');
    uiresume
    return
   case 'ok'
    % Go through the popup lists to make sure no marker is chosen
    % twice.

    
    markerids=get(figh,'UserData');
    poph=findobj(figh,'tag',markerids{1});
    mnames=get(poph,'string');
    
    mlist=cell(size(markerids));
    for i=1:length(markerids)
      poph=findobj(figh,'tag',markerids{i});
      mlist{i}=mnames{get(poph,'value')};
    end

    if (~(length(unique(mlist))==length(mlist)))
      set(stath,'ForegroundColor',[0.6 0 0]);
      set(stath,'String','Error: Markers multiply defined. Try again');
      return
    end
    
    set(figh,'UserData',cat(2,markerids,mlist));
    uiresume
  end
end


function notok=checkargs(v)
  notok=0;
  if (size(v{1},2)<2)
    %Check that the first argument is a n x 2 cell array.
    errordlg({'Wrong size of input argument ''markerid'' in',...
	      'function identify_markers.'})
    notok=1;
  elseif (size(v{1},1)>length(v{3}))
    % Check if set of markers in tsv is less than the markers
    % to be identified.
    errordlg({'Too few markers in tsv file.',...
	      'Error occured in function identify_markers.'})
    notok=1;
  end
  
  
function testfunction

  %test='left';
  test='right';
  
  if strcmp(test,'left')
    imname='eqdist_markers_left.png';
    [im,map]=imread(imname,'png');
    if (size(im,3)==1)
      im=ind2rgb(im,map);
    end
    
    markerid={'mc_medial',[0 0.98];...
	      'mc_lateral',[1 0.98];...
	      'mc_dorsal',[1 0.71];...
	      'fet_medial',[0 0.45];...
	      'fet_lateral',[1 0.45];...
	      'past_dorsal',[0 0.335];...
	      'past_medial',[0 0.25];...
	      'past_lateral',[1 0.31];...
	      'cof_medial',[0 0.16];...
	      'cof_lateral',[1 0.22];...
	      'cof_dorsal',[0 0.06];...
	      'hoof_lateral',[1 0.13];...
	      'hoof_dorsal',[1 0.05]};
    
    mnames={'mc_med_LH',...
	    'mc_lat_LH',...
	    'mc_dors_LH',...
	    'fet_med_LH',...
	    'fet_lat_LH',...
	    'past_med_LH',...
	    'past_lat_LH',...
	    'past_dors_LH',...
	    'cof_med_LH',...
	    'cof_lat_LH',...
	    'cof_dors_LH',...
	    'hoof_lat_LH',...
	    'hoof_dors_LH'};
  elseif strcmp(test,'right')
    imname='eqdist_markers_right.png';
    [im,map]=imread(imname,'png');
    if (size(im,3)==1)
      im=ind2rgb(im,map);
    end
    
    markerid={'mc_medial',[1 0.98];...
	      'mc_lateral',[0 0.98];...
	      'mc_dorsal',[0 0.71];...
	      'fet_medial',[1 0.45];...
	      'fet_lateral',[0 0.45];...
	      'past_dorsal',[1 0.335];...
	      'past_medial',[1 0.25];...
	      'past_lateral',[0 0.31];...
	      'cof_medial',[1 0.16];...
	      'cof_lateral',[0 0.22];...
	      'cof_dorsal',[1 0.06];...
	      'hoof_lateral',[0 0.13];...
	      'hoof_dorsal',[0 0.05]};
    
    mnames={'mc_med_RF',...
	    'mc_lat_RF',...
	    'mc_dors_RF',...
	    'fet_med_RF',...
	    'fet_lat_RF',...
	    'past_med_RF',...
	    'past_lat_RF',...
	    'past_dors_RF',...
	    'cof_med_RF',...
	    'cof_lat_RF',...
	    'cof_dors_RF',...
	    'hoof_lat_RF',...
	    'hoof_dors_RF'};
    
  end
  
  alist=identify_markers(markerid,im,mnames)
  
 
	  
	    
		    