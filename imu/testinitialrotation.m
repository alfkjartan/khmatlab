function testinitialrotation(pd, ntfrs, q12, p12);


  %% Check the initial rotation.
  q21 = qinv(q12);

  w2test = qrot(pd.prn.gyroCal', q21);
  a2test = qrot(pd.prn.accCal', q21);
  
  a2test2 = zeros(size(pd.pln.gyroCal))';
  alpha1 = centraldiff(pd.pln.gyroCal, 102.4);

  for i=1:size(pd.pln.gyroCal,1)
    a2test2(:,i) = slave_acc(pd.pln.accCal(i,:)', \
			     pd.pln.gyroCal(i,:)', alpha1(i,:)', p12);
  end

  figure(3)
  clf
  plot(pd.pln.timeCal(1:ntfrs), \
       pd.pln.accCal(1:ntfrs,1), 'color', [0 0 0.5])
  hold on
  plot(pd.pln.timeCal(1:ntfrs), \
	 pd.pln.accCal(1:ntfrs,2), 'color', [0 0.5 0])
    plot(pd.pln.timeCal(1:ntfrs), \
	 pd.pln.accCal(1:ntfrs,3), 'color', [0.5 0 0])
    plot(pd.prn.timeCal(1:ntfrs), \
	a2test(1,1:ntfrs)', 'b')
    plot(pd.prn.timeCal(1:ntfrs), \
	a2test(2,1:ntfrs)', 'g')
    plot(pd.prn.timeCal(1:ntfrs), \
	a2test(3,1:ntfrs)', 'r')
    plot(pd.pln.timeCal(1:ntfrs), \
	 a2test2(:, 1:ntfrs)', 'linewidth', 3);

    figure(5)
    clf
    plot(pd.pln.timeCal(1:ntfrs), \
	 pd.pln.gyroCal(1:ntfrs,1), 'color', [0 0 0.5])
    hold on
    plot(pd.pln.timeCal(1:ntfrs), \
	 pd.pln.gyroCal(1:ntfrs,2), 'color', [0 0.5 0])
    plot(pd.pln.timeCal(1:ntfrs), \
	 pd.pln.gyroCal(1:ntfrs,3), 'color', [0.5 0 0])
    plot(pd.prn.timeCal(1:ntfrs), \
	w2test(1,1:ntfrs)', 'b')
    plot(pd.prn.timeCal(1:ntfrs), \
	w2test(2,1:ntfrs)', 'g')
    plot(pd.prn.timeCal(1:ntfrs), \
	w2test(3,1:ntfrs)', 'r')
