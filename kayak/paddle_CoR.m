%% Script for testing how stationary the center of rotation is for a kayak paddle

%% Kjartan Halvorsen
%% 2013-05-15

plotit = 1;

n_parts = 6; %% The number of parts to divide each pull into.

skip_cycles = 3; %% Skip this number of cycles at beginning of trial

 
%% Data to use
dtapth = "/home/kjartan/Dropbox/Ergometer\ vs\ Wireless/mocap110416/tsvdata";

subjects = {"Albert", "Fredrik", "Joel"};

files = { "Fredrik/80_1.tsv",
	  "Fredrik/120_1.tsv",
	  "Fredrik/160_1.tsv",
	  "Joel/120_3.tsv"};


markers = {'L_shaft', 'L_sensor_front', 'L_sensor_back', 'R_shaft', 'R_sensor_front', 'R_sensor_back'};

%for i=1:length(files)
for i=3:3
    f = files{i};
    md = openmocapfile('', fullfile(dtapth, f));

    
    pd = extractmarkers(md, markers); % Nframes x Nmarkers*3
    
    %%pdd = zeros(3, size(pd,1), length(markers));
    %%for i=1:length(markers)
    %%pdd(:,:,i) = pd(:,(i-1)*3+1:i*3)';
    %%end

    %% Find start of pull phase on left and right side 
    l_shaft = extractmarkers(md, {'L_shaft'});
    r_shaft = extractmarkers(md, {'r_shaft'});

    freq = str2num(getvalue(md{1}, 'FREQUENCY'));



    [left_starts, left_ends] = peakdet(l_shaft(:,1), 0.03);
    %% Skip the first starts, and any end before first start
    left_starts = left_starts(find(left_starts>left_starts(skip_cycles)));
    left_ends = left_ends(find(left_ends>left_starts(1)));

    cycles_l = min(length(left_starts), length(left_ends));
    left_starts = left_starts(1:cycles_l);
    left_ends = left_ends(1:cycles_l);
    
    [right_starts, right_ends] = peakdet(r_shaft(:,1), 0.03);
    %% Skip the first starts, and any end before first start
    right_starts = right_starts(find(right_starts>right_starts(skip_cycles)));
    right_ends = right_ends(find(right_ends>right_starts(1)));

    cycles_r = min(length(right_starts), length(right_ends));
    right_starts = right_starts(1:cycles_r);
    right_ends = right_ends(1:cycles_r);
    
    %% Limit the number of cycles during coding
    %%left_starts = left_starts(1:8);
    %%left_ends = left_ends(1:8);

    % There are about 100 frames in each pull. Divide into n_parts parts
    left_CoR = zeros(3,n_parts, length(left_starts));
    left_force = zeros(n_parts, length(left_starts));
    for c = 1:length(left_starts)
	%% Compute velocity at end of shaft
	inds = left_starts(c):left_ends(c);
	dist2flywheel = sqrt(sum((l_shaft(inds, :) - repmat(l_shaft(inds(1),:), length(inds), 1)).^2, 2));
	paddlevel = centraldiff(dist2flywheel, freq);
	pull_length = (left_ends(c) - left_starts(c))*0.8;
	part_length = floor(pull_length/n_parts);
	for p = 1:n_parts
	  %%cycledata = pdd(:, (left_starts(c)+(p-1)*part_length):(left_starts(c)+p*part_length),:);
	    p_inds = (left_starts(c)+(p-1)*part_length):(left_starts(c)+p*part_length);
	    v_inds = p_inds - inds(1) + 1;
	    cycledata = pd(p_inds, :);
	    left_CoR(:, p, c) = ehrig_jc(cycledata);
	    %% Force is approximately proportional to acceleration of flywheel
	    left_force(p, c) = ( mean(paddlevel(v_inds(end-3:end))) - mean(paddlevel(v_inds(1:3))) ) / length(v_inds) * freq;
	end
    end
    %% Normalize force 
    left_force = left_force / max(max(abs(left_force)));


    right_CoR = zeros(3,n_parts, length(right_starts));
    right_force = zeros(n_parts, length(right_starts));
    for c = 1:length(right_starts)
	%% Compute velocity at end of shaft
	inds = right_starts(c):right_ends(c);
	dist2flywheel = sqrt(sum((r_shaft(inds, :) - repmat(r_shaft(inds(1),:), length(inds), 1)).^2, 2));
	paddlevel = centraldiff(dist2flywheel, freq);
	pull_length = (right_ends(c) - right_starts(c))*0.8;
	part_length = floor(pull_length/n_parts);
	for p = 1:n_parts
	    p_inds = (right_starts(c)+(p-1)*part_length):(right_starts(c)+p*part_length);
	    v_inds = p_inds - inds(1) + 1;
	    cycledata = pd(p_inds,:);
	    right_CoR(:, p, c) = ehrig_jc(cycledata);
	    %% Force is approximately proportional to acceleration of flywheel
	    right_force(p, c) = ( mean(paddlevel(v_inds(end-3:end))) - mean(paddlevel(v_inds(1:3))) ) / length(v_inds) * freq;
	end
    end
    %% Normalize force 
    right_force = right_force / max(max(abs(right_force)));

    STD = 1;     %% Use this number of standard deviations
    conf = 2*normcdf(STD) - 1;  %% Covers about 95% of the population
    scale = chi2inv(conf, 3);  %% Inverse chi squared
    mu_left = mean(reshape(left_CoR, 3, length(left_starts)*n_parts)')';
    cov_left = scale * cov(reshape(left_CoR, 3, length(left_starts)*n_parts)');
    [V_left,D_left] = eig(cov_left);
    [D_left, order] = sort(diag(D_left), 'descend');
    D_left = diag(D_left);
    V_left = V_left(:, order);

    mu_right = mean(reshape(right_CoR, 3, length(right_starts)*n_parts)')';
    cov_right = scale * cov(reshape(right_CoR, 3, length(right_starts)*n_parts)');
    [V_right,D_right] = eig(cov_right);
    [D_right, order] = sort(diag(D_right), 'descend');
    D_right = diag(D_right);
    V_right = V_right(:, order);
    if plotit
      figure(1)
      clf
      plot(l_shaft(:,1))
      hold on

      yl = get(gca, 'ylim');
      for s = left_starts
	  plot([s s], yl, 'g')
      end
      for s = left_ends
	  plot([s s], yl, 'r')
      end

      figure(4)
      clf
      plot(r_shaft(:,1))
      hold on

      yl = get(gca, 'ylim');
      for s = right_starts
	  plot([s s], yl, 'g')
      end
      for s = right_ends
	  plot([s s], yl, 'r')
      end

      figure(2)
      clf
      for c = 1:length(left_starts)
	inds = left_starts(c):left_ends(c);
	plot(l_shaft(inds,1), l_shaft(inds,2), 'r')
	hold on
	for p = 1:n_parts
	    plot(left_CoR(1, p, c), left_CoR(2, p, c), 'ro', 'linewidth', 1, 'markersize', 14*left_force(p,c))
	end
      end
      for c = 1:length(right_starts)
	inds = right_starts(c):right_ends(c);
	plot(r_shaft(inds,1), r_shaft(inds,2), 'g')
	for p = 1:n_parts
	    plot(right_CoR(1, p, c), right_CoR(2, p, c), 'go', 'linewidth', 1, 'markersize', 14*right_force(p,c))
	end
      end

      t = linspace(0, 2*pi, 200);
      ee = [cos(t); sin(t)];  % Unit circle
      VV_left = V_left*sqrt(D_left);
      VV_right = V_right*sqrt(D_right);
     
      eee = VV_left(1:2, 1:2)*ee + repmat(mu_left(1:2), 1, 200);
      plot(eee(1,:), eee(2,:), 'color', [0.5 0 0], 'linewidth', 2)
      eee = VV_right(1:2, 1:2)*ee + repmat(mu_right(1:2), 1, 200);
      plot(eee(1,:), eee(2,:), 'color', [0 0.5 0], 'linewidth', 2)
      xlabel('x [m]')
      ylabel('y [m]')
      title(f)
      axis equal
      fname = join(split(f, '/'), '_')
      print([fname, '_horizontal.pdf'], '-dpdf')


      figure(3)
      clf
      for c = 1:length(left_starts)
	inds = left_starts(c):left_ends(c);
	plot(l_shaft(inds,1), l_shaft(inds,3), 'r')
	hold on
	for p = 1:n_parts
	    plot(left_CoR(1, p, c), left_CoR(3, p, c), 'ro', 'linewidth', 1, 'markersize', 14*left_force(p,c))
	end
      end
      for c = 1:length(right_starts)
	inds = right_starts(c):right_ends(c);
	plot(r_shaft(inds,1), r_shaft(inds,3), 'g')
	for p = 1:n_parts
	    plot(right_CoR(1, p, c), right_CoR(3, p, c), 'go', 'linewidth', 1, 'markersize', 14*right_force(p,c))
	end
      end

      eee = VV_left([1 3], [1 3])*ee + repmat(mu_left([1 3]), 1, 200);
      plot(eee(1,:), eee(2,:), 'color', [0.5 0 0], 'linewidth', 2)
      eee = VV_right([1 3], [1 3])*ee + repmat(mu_right([1 3]), 1, 200);
      plot(eee(1,:), eee(2,:), 'color', [0.0 0.5 0.], 'linewidth', 2)


      xlabel('x [m]')
      ylabel('z [m]')
      title(f)
      axis equal
      fname = join(split(f, '/'), '_')
      print([fname, '_saggital.pdf'], '-dpdf')

    end
end

    
