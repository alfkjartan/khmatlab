function get_hand_force(ref, pliancepos, mpdata, pointerdist)
  %% get_hand_force(ref, pliancepos, mdata, pdata)
  %% Computes the resulting grip force on the handle, given the pliance
  %% data and marker data for the grip.
  %%
  %% Input
  %%   ref         ->  Reference marker  positions
  %%   pliancepos  ->  Cell array with data for the position of the
  %%                   pliance sensor:
  %%                   {[x,y], mdatafilename; ...}
  %%   mpdata      ->  Cell array with file names, ntrials x 2  with
  %%                   {markerdatafile1, pliancedatafile1
  %%                    markerdatafile2, pliancedatafile2,...}
  %%   pdata       ->  Cell array with file names for pliance data

  %% Kjartan Halvorsen
  %% 2012-02-09

  %% Solution:
  %% 1) The marker data in the pliancepos cell array gives the position, in
  %% 3d, of some of the cells in the pliance mat. 
  %% 2) A cylinder of known diameter (20mm) with known positions of all
  %% the cells on the mat is fitted to the data, based on the knowledge
  %% of the position of some of the cells. This is done simply by using
  %% soder.m (point correspondance).
  %% 3) The model of the pliance cylinder contains a  list of all the
  %% cells, their position in a local coordinate system, and their
  %% normal vector. The model is loaded with get_pliance_model.m

  %% virtual markers at all the cells.

  cylindernames = {'cylinder_tip';'cylinder_left';'cylinder_right'};
  pointernames = {'pointer_tip'; 'pointer_handle'};

  % Find reference positions for the cylinder markers
  refd = openmocapfile('', ref);
  cylref = extractmeanmarkers(refd, cylindernames);

  % Go through the set of pliance data and find positions for the cells
  npoints = length(pliancepos);
  cellposdata = nan(3,npoints);
  cellpos = nan(npoints,1);
  for i=1:npoints
    pd = openmocapfile('', pliancepos{i});

    pd{2} = pd{2}(100:150,:); % Use only part of the data.
    
    [dr,name] = fileparts(pliancepos{i}); 
    [tokens] = regexp(name, '^pos(3|9|13)([1-4])$', 'tokens');
    %%keyboard
    %%row = str2num(tokens{1}{2});
    %%col = str2num(tokens{1}{3});
    row = str2num(tokens{1}{1});
    col = str2num(tokens{1}{2});
    cellpos(i) = (row-1)*4 + col;

    pointerd = extractmeanmarkers(pd, pointernames);
    v = pointerd(:,1) - pointerd(:,2);
    v = v / norm(v);
    tip = pointerd(:,1) + v*pointerdist;
    cpd = [tip;1];
    
    cdt = extractmeanmarkers(pd, cylindernames);
    
    T = soder(cat(1, reshape(cdt, 1, 9), reshape(cylref, 1, 9)));
    tipref = T*cpd;
    cellposdata(:,i) = tipref(1:3);
  end

  %% Get the pliance grip model, and determine the transformation that
  %% takes points in the local coordinate system of the pliance handle to
  %% the ref position.
  [plpoints, plnvectors, plnames, cellarea] = get_pliance_model();
  %%keyboard
  plps = plpoints(:,cellpos);
  Tp = soder(cat(1, reshape(plps(:), 1, npoints*3), ...
		 reshape(cellposdata(:), 1, npoints*3)));
  gripref = Tp*cat(1, plpoints, ones(1, size(plpoints,2)));
  gripnormalref = Tp*cat(1, plnvectors, zeros(1, size(plnvectors,2)));


  %%  Go through the trials. Create virtual points for all the sensors and
  %%  save to tsv files.
  %%  Load corresponding pliance file 
  
  for i=1:size(mpdata, 1)
    md = openmocapfile('', mpdata{i,1});
    nframes = size(md{2}, 1);
    T = getMotion(refd, md, cylindernames);
    %grippoints = movePoints(T, gripref);
    %gripnvecs = movePoints(T, gripnormalref);

    %% Get grip force and cop in local frame of grip
    [gripf_loc, gripcop_loc] = get_grip_force(mpdata{i,2});

    figure(10)
    clf
    plot(gripf_loc')
    title('Local grip force [N]')
    hold on
    gripf_loc_magn = sqrt(sum(gripf_loc.^2));
    plot(gripf_loc_magn, 'linewidth', 3)
    ylabel('N')
    xlabel('Frame')
    print -depsc2 force_local_ref.eps
    system('epstopdf force_local_ref.eps')

    figure(11)
    clf
    plot(gripcop_loc')
    title('Local grip CoP [m]')
    print -depsc2 CoP_local_ref.eps
    system('epstopdf CoP_local_ref.eps')

    %%keyboard
    
    nframes = min(size(gripf_loc, 2), size(md{2}, 1));
    
    gripforce = zeros(3, nframes);
    gripcop = zeros(3,nframes);
    grippoints = zeros(3, 64, nframes);
    gripnvecs = zeros(3, 64, nframes);

    for j = 1:nframes
      Tj = reshape(T(j,:), 4, 4);
      gripf = Tj * Tp*[gripf_loc(:,j); 0];
      gripforce(:,j) = gripf(1:3);
      gripc = Tj * Tp* [gripcop_loc(:,j); 1];
      gripcop(:,j) = gripc(1:3);
      gripp = Tj*gripref;
      grippoints(:,:,j) = gripp(1:3,:);
      gripn = Tj*gripnormalref;
      gripnvecs(:,:,j) = gripn(1:3,:);
    end

    figure(1)
    clf
    plot(gripforce')
    hold on
    plot(sqrt(sum(gripforce.^2)), 'linewidth', 3)
    title('Force')
    ylabel('N')
    xlabel('Frame')
    print -depsc2 force_global_ref.eps
    system('epstopdf force_global_ref.eps')
    figure(2)
    clf
    plot(gripcop')
    title('CoP')

    figure(4)
    clf
    plot(T)
    
    figure(5)
    clf
    plot(gripf_loc')
    hold on
    plot(sqrt(sum(gripf_loc.^2)), 'linewidth', 3)
    title('Local grip force')
    
    vnames = plnames;
    for j = 1:length(plnames)
      vnames{j} = [vnames{j}, "_vec"];
    end

    % Write tsv
    [pth, name, ext] = fileparts(mpdata{i,1});
    fnameout = fullfile(pth, [name, '_grip', ext]);
    mnames = getvalue(md{1}, 'MARKER_NAMES');
    mnames = cat(1, mnames, {'grip_force'; 'grip_CoP'}, plnames, vnames);
    fid = fopen(fnameout, 'w');
    md{1} = putvalue(md{1}, 'NO_OF_MARKERS', sprintf("%d", length(mnames)));
    grippointstsv = reshape(grippoints(:), length(plnames)*3, nframes)'; 
    gripnvecstsv = grippointstsv + 0.01*reshape(gripnvecs(:), length(plnames)*3, nframes)'; 
    write3dtsv(putvalue(md{1}, 'MARKER_NAMES', mnames),...
	       cat(2, md{2}, gripcop'+gripforce'/10, gripcop',...
		   grippointstsv, gripnvecstsv)*1000, fid);
    fclose(fid);
    disp(['Wrote file ', fnameout])
    keyboard
  end

  
  
  
  
