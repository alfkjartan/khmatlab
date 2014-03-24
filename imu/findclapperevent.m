function [fr0_mc, fr0_prn, fr0_pln] = findclapperevent(clap1, clap2, \
						       accR, accL)
%% Returns the fram index of the same clapper event. Found by visual inspection
    
  v = clap2 - clap1;
  vclap2 = centraldiff(clap2, 1);
  vv = sum(v.*vclap2,2);

  figure(1)
  clf
  plot(vv)

  %%accL = sum(accL.^2, 2);
  %%accR = sum(accR.^2, 2);
  %%accL = accL(:,3);
  %%accR = accR(:,2);

  figure(2)
  clf
  plot(accL)


  figure(3)
  clf
  plot(accR)

  keyboard