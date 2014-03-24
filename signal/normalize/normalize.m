function normalize
  % Program that normalizes a time series to 0 to 100
  
  % Kjartan Halvorsen
  % 2004-02-27

  oldd = pwd;
  
  preandpostdata = 11; % In procent

  answer = 'Yes';

  figh=1;
  
  while (strcmp(answer, 'Yes'))
    % Load data
    [mdata, slask, filename] = opentsv('Open the data file');

    if ~isempty(mdata)
  
      nfrs = size(mdata{2},1);
      eventfr = fix([0.2*nfrs; 0.8*nfrs]);
      [eventfr,figh]=detectstepsplot(eventfr, mdata{2}, ...
				     'Mark start and stop',...
				     ['Adjust lines to mark start and stop',...
		    'of region to normalize']);
      close(figh);

      lngth = eventfr(end) - eventfr(1) + 1;
      percentlength = lngth / 100; % Lenght of 1 percent in samples

      startfr1 = eventfr(1) - fix(preandpostdata*percentlength);
      startfr = max(1, startfr1);
      startpercent = ceil( (startfr - eventfr(1))/ percentlength );
      startfr = floor( max(1, eventfr(1) + startpercent*percentlength));
      
      stopfr1 = eventfr(end) + fix(preandpostdata*percentlength);
      stopfr = min(nfrs, stopfr1);
      endpercent = 100 + floor(stopfr - eventfr(end))/percentlength;
      stopfr = min( nfrs, ...
		    eventfr(end) + ceil( (endpercent-100)*percentlength ));

      timelinenew = (startpercent:endpercent);
      timelineoriginal = ((startfr:stopfr) - eventfr(1)) / percentlength;
      
      ndata = interp1(timelineoriginal', mdata{2}(startfr:stopfr,:),...
		      timelinenew);
      
      figh=plot(timelinenew',ndata);
      
      
      % Save data

      [pthpart, namepart, ext] = fileparts(filename);
      
      cd(pthpart);
      
      fid = fopen( fullfile(pthpart, [namepart, '_norm.txt']),'w' );

      if (~isempty(mdata{1}))
	header = (mdata{1})';
	header1 = header(:,1:end-1);
	header2 = header(:,end);
	if ~isempty(header1)
	  fprintf(fid, '%s\t', header1{:});
	fprintf(fid,'\n');
	end
	fprintf(fid, '%s\t', 'time', header2{:});
	fprintf(fid,'\n');
      end
      ok=write3dtsv(cat(2, timelinenew', ndata), fid);
      fclose(fid);
    end
    
    answer = questdlg('Normalize another data file?',...
		      'Question', 'Yes', 'No', 'Yes');
    
  end
  
  
  cd(oldd);
  
  
  
  