function [command,item] = saveobj(object)
%SAVEOBJ stores the information necessary to recreate an object.
% The object can be a line, patch, surface, text, or image.
% Neither axes nor figures can be stored using SAVEOBJ.
% [COMMAND,DATA] = SAVEOBJ(OBJECT)
% returns a COMMAND string and a structure DATA which 
% keeps information about object properties that cannot be stored in
% strings. The object may then be recreated using the
% LOADOBJ command. 

% Copyright (c) 1997 by Keith Rogers


if(~ishandle(object))
	error('Not an object!');
end
UserData = get(object,'UserData');
Type = get(object,'Type');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following are common to all objects and %
% can be fully specified by text strings.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CommonBlock = ['''ButtonDownFcn'',''' get(object,'ButtonDownFcn') ''',' ...
               '''CreateFcn'',''' get(object,'CreateFcn') ''',' ...
               '''DeleteFcn'',''' get(object,'DeleteFcn') ''',' ...
               '''HandleVisibility'',''' get(object,'HandleVisibility') ''',' ...
			   '''Clipping'',''' get(object,'Clipping') ''',' ...
			   '''Interruptible'',''' get(object,'Interruptible') ''',' ...
			   '''SelectionHighlight'',''' get(object,'SelectionHighlight') ''',' ...
			   '''Tag'',''' get(object,'Tag') ''',' ...
			   '''Selected'',''' get(object,'Selected') ''',' ...
			   '''Visible'',''' get(object,'Visible') ''','];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff unique to LINE objects %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(Type)
	case 'line',
		command = ['line(' ...
				   '''Color'',Color,' ...
				   CommonBlock ...
				   '''UserData'',UserData);'];
		item(1).val = UserData;
		item(1).prop = 'UserData';
		item(2).val = get(object,'Color');
		item(2).prop = 'Color';

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Stuff unique to SURFACE objects %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'surface',
		
		command = ['surface(' ...
				   '''MeshStyle'',''' get(object,'MeshStyle') ''',' ...			   
				   CommonBlock ...
				   '''UserData'',UserData);'];
		item(1).val = UserData;
		item(1).prop = 'UserData';

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Stuff unique to TEXT objects %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'text',
		command = ['text(' ...
				   '''Color'',Color,' ...
				   '''EraseMode'',''' get(object,'EraseMode') ''',' ...
				   '''FontAngle'',''' get(object,'FontAngle') ''',' ...
				   '''FontName'',''' get(object,'FontName') ''',' ...
				   '''FontSize'',' num2str(get(object,'FontSize')) ',' ...			   			   
				   '''FontWeight'',''' get(object,'FontWeight') ''',' ...
				   '''HorizontalAlignment'',''' get(object,'HorizontalAlignment') ''',' ...
				   '''Position'',Position' ',' ...			   			   
				   '''Rotation'',' num2str(get(object,'Rotation')) ',' ...			   			   
				   '''String'',''' get(object,'String') ''',' ...
				   '''Units'',''' get(object,'Units') ''',' ...
				   '''VerticalAlignment'',''' get(object,'VerticalAlignment') ''',' ...
				   CommonBlock ...
				   '''UserData'',UserData);'];
		item(1).val = UserData;
		item(1).prop = 'UserData';
		item(2).val = get(object,'Color');
		item(2).prop = 'Color';
		item(3).val = get(object,'Position');
		item(3).prop = 'Position';

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Stuff unique to PATCH objects %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'patch',
		command = ['patch(' ...
				   '''Vertices'',Vertices,' ...			   			   
				   '''FaceVertexCData'',FaceVertexCData' ',' ...			   			   
				   '''Faces'',Faces,' ...
				   CommonBlock ...
				   '''UserData'',UserData);'];
		item(1).val = UserData;
		item(1).prop = 'UserData';
		item(2).val = get(object,'Vertices');
		item(2).prop = 'Vertices';
		item(3).val = get(object,'FaceVertexCData');
		item(3).prop = 'FaceVertexCData';
		item(4).val = get(object,'Faces');
		item(4).prop = 'Faces';

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Stuff unique to IMAGE objects %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'image',
		command = ['image(' ...
				   '''XData'',XData,' ...			   			   
				   '''YData'',YData,' ...			   			   
				   '''CData'',CData,' ...
				   '''CDataMapping'',''' get(object,'CDataMapping') ''',' ...			   
				   '''CDataScaling'',''' get(object,'CDataScaling') ''',' ...			   
				   '''EraseMode'',''' get(object,'EraseMode') ''',' ...			   
				   CommonBlock ...
				   '''UserData'',UserData);'];
		item(1).val = UserData;
		item(1).prop = 'UserData';
		item(2).val = get(object,'XData');
		item(2).prop = 'XData';
		item(3).val = get(object,'YData');
		item(3).prop = 'YData';
		item(4).val = get(object,'CData');
		item(4).prop = 'CData';
end			   	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff common to lines, patches,  and surfaces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(Type,'line')|strcmp(Type,'surface')|strcmp(Type,'patch'))
	command(end-1:end) = ', ';
	command = [command ...
	   '''EraseMode'',''' get(object,'EraseMode') ''',' ...
	   '''LineStyle'',''' get(object,'LineStyle') ''',' ...			   
	   '''Marker'',''' get(object,'Marker') ''',' ...			   
	   '''Linewidth'',' num2str(get(object,'LineWidth')) ',' ...			   			   
	   '''MarkerSize'',' num2str(get(object,'MarkerSize')) ',' ...			   			   
	   '''MarkerEdgeColor'',MarkerEdgeColor' ',' ...			   			   
	   '''MarkerFaceColor'',MarkerFaceColor' ',' ...			   			   
	   '''XData'',XData' ',' ...			   			   
	   '''YData'',YData' ',' ...			   			   
	   '''ZData'',ZData' ');']; 
	item(end+1).val = get(object,'MarkerEdgeColor');
	item(end).prop = 'MarkerEdgeColor';
	item(end+1).val = get(object,'MarkerFaceColor');
	item(end).prop = 'MarkerFaceColor';
	item(end+1).val = get(object,'XData');
	item(end).prop = 'XData';
	item(end+1).val = get(object,'YData');
	item(end).prop = 'YData';
	item(end+1).val = get(object,'ZData');
	item(end).prop = 'ZData';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff common to patches and surfaces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(Type,'patch')|strcmp(Type,'surface'))
	command(end-1:end) = ', ';
	command = [command ...
	   '''EdgeColor'',EdgeColor,' ...
	   '''FaceColor'',FaceColor,' ...
	   '''CData'',CData,' ...
	   '''CDataMapping'',''' get(object,'CDataMapping') ''',' ...
	   '''FaceLighting'',''' get(object,'FaceLighting') ''',' ...
	   '''EdgeLighting'',''' get(object,'EdgeLighting') ''',' ...			   
	   '''BackFaceLighting'',''' get(object,'BackFaceLighting') ''',' ...			   
	   '''NormalMode'',''' get(object,'NormalMode') ''',' ...			   
	   '''AmbientStrength'',' num2str(get(object,'AmbientStrength')) ',' ...			   			   
	   '''DiffuseStrength'',' num2str(get(object,'DiffuseStrength')) ',' ...			   			   
	   '''SpecularStrength'',' num2str(get(object,'SpecularStrength')) ',' ...			   			   
	   '''SpecularExponent'',' num2str(get(object,'SpecularExponent')) ',' ...			   			   
	   '''SpecularColorReflectance'',' num2str(get(object,'SpecularColorReflectance')) ',' ...			   			   
	   '''VertexNormals'',VertexNormals' ');']; 
	item(end+1).val = get(object,'EdgeColor');
	item(end).prop = 'EdgeColor';
	item(end+1).val = get(object,'FaceColor');
	item(end).prop = 'FaceColor';
	item(end+1).val = get(object,'VertexNormals');
	item(end).prop = 'VertexNormals';
	item(end+1).val = get(object,'CData');
	item(end).prop = 'CData';
end
