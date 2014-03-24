function [pd, md, time_mc, fr0_pln, fr0_prn, fr0_mc] = get_calibrated_validationdata(fn_pln, fn_prn, fn_mc, ...
						  fn_sync);

%% Will load data and return calibrated and time-shifted data. If
%% indices of starting frames are not given, then this must be done by
%% visual inspection. Plots are generated for this.

  plotit = 0;

  synchronize = (nargin < 4);


  pd.prn = opensensordata(fn_prn);
  pd.pln = opensensordata(fn_pln);

  %% prn and pln in turn contains the fields
  %%   raw
  %%   time
  %%   forceBin
  %%   accBin
  %%   gyroBin

  md = openmocapfile('', fn_mc);

  sfreq =  str2double(getvalue(md{1}, 'FREQUENCY'));

  
  if (synchronize)
    clap1 = extractmarkers(md, {'clap1'});
    clap2 = extractmarkers(md, {'clap2'});
    clap1(find(clap1==0)) = NaN;
    clap2(find(clap2==0)) = NaN;
    findclapperevent(clap1, clap2, \
		     pd.prn.accBin(:,2), \
		     pd.pln.accBin(:,3));
    return
  else
    [fr0_pln, fr0_prn, fr0_mc] = feval(fn_sync);
  end

  %% Use found event in PLN
  plntime = pd.pln.time(fr0_pln);
  prnbefore = find(pd.prn.time < plntime);
  %% Take the next frame of data in PRN 
  fr0_prn = prnbefore(end) + 1;

  %% sync between nodes:
  disp('Time pln at clap:')
  disp(pd.pln.time(fr0_pln))
  
  disp('Time prn at clap:')
  disp(pd.prn.time(fr0_prn))
  
  if plotit
    %% Look at the data
    nfr = 1000;
    
    figure(3)
    clf
    t = (-50:nfr)*0.01;
    plot(t, pd.pln.gyroBin(fr0_pln-50:fr0_pln+nfr,:))

    figure(4)
    clf
    plot(t, pd.prn.gyroBin(fr0_prn-50:fr0_prn+nfr,:))
  end  
  
  %% Calibration data. Linear model for all sensors:
  %% SI_unit = cal*(binary-offset);
  
  %% Load calibration data
  acc_calibration; %% 

  %% From Freescale accel documentation: (Not used)
  %%mV_per_g = 180;
  %%mV_per_bit = 3.3*1000/1024;
  %%g_per_bit = 1/mV_per_g*mV_per_bit;
  %%g_per_bit = 1/mV_per_g*mV_per_bit*0.6; % Including fudge factor
  apbL = 9.825./pln_acc_bits_per_g;
  apbR = 9.825./prn_acc_bits_per_g;

  acc0L = pln_acc_zero';
  accCalL = diag(apbL); %% m/s/s/bit
  
  acc0R = prn_acc_zero';
  accCalR = diag(apbR); %% m/s/s/bit
  

  %% Gyro calibration
  gyro_calibration;
  
  %% Invensense IDG-500 gyro. DATA NOT USED
  %%mV_per_deg_per_s = 2.0*3.3/3.0;
  %%mV_per_bit = 3.3*1000/1024;
  %%rad_per_s_per_bit = (pi/180.0) * 1.0/mV_per_deg_per_s * mV_per_bit;
  %%rpspb = rad_per_s_per_bit;
  gyroCalL = -diag(1.0 ./ pln_raw_per_rad_per_s); %% OBS: negative to
  %% switch direction of the gyro axes.
  gyro0L = pln_gyro_zero';
  
  gyroCalR = -diag(1.0 ./ prn_raw_per_rad_per_s);
  gyro0R = prn_gyro_zero';
  
  %% Apply calibration 
  timeCalL = ( pd.pln.time(fr0_pln:end) - pd.pln.time(fr0_pln) ) ; % *1.024;
  nfrLc = timeCalL(end)*100;
  nfrL = length(timeCalL);
  pd.pln.timeCal = timeCalL/1000;
  
  pd.pln.accCal = (pd.pln.accBin(fr0_pln:end,:) ...
      - repmat(acc0L, nfrL, 1) ) * accCalL;
  
  pd.pln.gyroCal = (pd.pln.gyroBin(fr0_pln:end,:) ...
     - repmat(gyro0L, nfrL, 1) ) * gyroCalL;

  %%pd.pln.accCal = interp1(timeCalL, ( pd.pln.accBin(fr0_pln:end,:) - \
  %%				     repmat(acc0L, nfrL, 1) ) * accCalL, pd.pln.timeCal);
  
  %%pd.pln.gyroCal = interp1(timeCalL, ( pd.pln.gyroBin(fr0_pln:end,:) - \
  %%				      repmat(gyro0L, nfrL, 1) ) * gyroCalL, pd.pln.timeCal);
  %%
  timeCalR = ( pd.prn.time(fr0_prn:end) - pd.prn.time(fr0_prn) );
  nfrRc = timeCalR(end)*100;
  nfrR = length(timeCalR);

  %%pd.prn.timeCal = linspace(timeCalR(1), timeCalR(end), nfrRc+1);
  pd.prn.timeCal = timeCalR/1000;
  pd.prn.accCal = (pd.prn.accBin(fr0_prn:end,:) ...
      - repmat(acc0L, nfrR, 1) ) * accCalR;
  pd.prn.gyroCal = (pd.prn.gyroBin(fr0_prn:end,:) ...
      - repmat(gyro0R, nfrR, 1) ) * gyroCalR;

 

  md{2} = md{2}(fr0_mc:end,:);
  nfr_mc = size(md{2},1);
  
  time_mc = (0:nfr_mc-1)/sfreq;
