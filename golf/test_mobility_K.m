%% Script for testing how the parametrization of the internal motion, K, influences 
%% mobility calculations. Should be independent!

%% Kjartan Halvorsen
%% 2013-11-26

%% Assumes    [gmleft, gmright, gmbothleft, gmbothright] = build_test_models(...)

test_double_hand = 1;
test_single_hand = 0;

NN = 10;

if test_double_hand
  gm = gmbothleft;
  nst = size(gm.gcnames,1);
  states = zeros(nst,2);

  Jh_ep = end_point_jacobian(gm, states);
  Mh = generalized_manipulator_inertia(gm, states);

  Gtransp_ep = repmat(eye(3), [2,1]);
  Gbartransp_ep = [Gtransp_ep  zeros(6, nst-2*3)
		   zeros(nst-2*3,3) eye(nst-2*3)];

  Kh = null(Jh_ep(:,:,1));
  Kh2 = 2*Kh;

  Jbar_ep = [Jh_ep(:,:,1); Kh'];
  Jbar2_ep = [Jh_ep(:,:,1); Kh2'];
  Jbarinv = inv(Jbar_ep);
  Jbarinv2 = inv(Jbar2_ep);

  Mbar1 = Gbartransp_ep'*Jbarinv'*Mh(:,:,1)*Jbarinv*Gbartransp_ep;

  Mbar2 = Gbartransp_ep'*Jbarinv2'*Mh(:,:,1)*Jbarinv2*Gbartransp_ep;

  tic();
  for i=1:NN
    Jbar_ep = [Jh_ep(:,:,1); Kh'];
    Jbar2_ep = [Jh_ep(:,:,1); Kh2'];
    Jbarinv = inv(Jbar_ep);
    Mbar1 = Gbartransp_ep'*Jbarinv'*Mh(:,:,1)*Jbarinv*Gbartransp_ep;
    W1 = inv(Mbar1);
  end
  toc()
  W2 = inv(Mbar2);

  assert(W1(1:3,1:3), W2(1:3,1:3), 1e-12)



  %% New way of computing mobility from paper (eq 24)
  K0 = Kh';
  JJ = Jh_ep(:,:,1);
  M = Mh(:,:,1);

  tic();
  for i=1:NN

    Jdagger = JJ' * inv(JJ*JJ');
    JdaggerGtransp = Jdagger*Gtransp_ep;

    F = JdaggerGtransp' * ( M - M * K0' * inv(K0*M*K0') * K0 * M) * JdaggerGtransp;

    W3 = inv(F);
  end
  toc()

  assert(W1(1:3, 1:3), W3, 1e-12)


  [W4slask, W4] = golfer_mobility(gm, states);

  assert(W3(1:3, 1:3), W4(1:3, 1:3), 1e-12)


  %% Test hypothesis that W = G*Jh*Yh*Jh'*G'
  %% FALSE!!
  W5 = Gtransp_ep'*JJ*inv(M)*JJ'*Gtransp_ep;
  assert(W3(1:3, 1:3), W5(1:3, 1:3), 1e-12)
  
  %% Test inv( M - M*K'*inv(K*M*K')*K'*M ) = inv(M)
  Minv2 = inv(M - M*K0'*inv(K0*M*K0')*K0*M);
  assert(inv(M), Minv2, 1e-12)


end

if test_single_hand
  %% Test that result is correct for one handed manipulation

  gm = gmleft;
  nst = size(gm.gcnames,1);
  states = zeros(nst,2);

  Jh_ep = end_point_jacobian(gm, states);
  Mh = generalized_manipulator_inertia(gm, states);

  Gtransp_ep = eye(3);
  Gbartransp_ep = [eye(3)  zeros(3, nst-3)
		   zeros(nst-3,3) eye(nst-3)];

  Kh = null(Jh_ep(:,:,1));
  Kh2 = 2*Kh;

  Jbar_ep = [Jh_ep(:,:,1); Kh'];
  Jbar2_ep = [Jh_ep(:,:,1); Kh2'];
  Jbarinv = inv(Jbar_ep);
  Jbarinv2 = inv(Jbar2_ep);


  Mbar1 = Gbartransp_ep'*Jbarinv'*Mh(:,:,1)*Jbarinv*Gbartransp_ep;

  Mbar2 = Gbartransp_ep'*Jbarinv2'*Mh(:,:,1)*Jbarinv2*Gbartransp_ep;

  tic();
  for i=1:NN
    Jbar_ep = [Jh_ep(:,:,1); Kh'];
    Jbar2_ep = [Jh_ep(:,:,1); Kh2'];
    Jbarinv = inv(Jbar_ep);
    Mbar1 = Gbartransp_ep'*Jbarinv'*Mh(:,:,1)*Jbarinv*Gbartransp_ep;
    W1 = inv(Mbar1);
  end
  toc()
  W2 = inv(Mbar2);

  assert(W1(1:3,1:3), W2(1:3,1:3), 1e-12)



  %% New way of computing mobility from paper (eq 24)
  K0 = Kh';
  JJ = Jh_ep(:,:,1);
  M = Mh(:,:,1);

  tic();
  for i=1:NN

    Jdagger = JJ' * inv(JJ*JJ');
    JdaggerGtransp = Jdagger*Gtransp_ep;

    F = JdaggerGtransp' * ( M - M * K0' * inv(K0*M*K0') * K0 * M) * JdaggerGtransp;

    W3 = inv(F);
  end
  toc()

  assert(W1(1:3, 1:3), W3, 1e-12)


  [W4slask, W4] = golfer_mobility(gm, states);

  assert(W3(1:3, 1:3), W4(1:3, 1:3), 1e-12)

  % Original way of computing
  W5 = JJ*inv(M)*JJ';
  assert(W4(1:3, 1:3), W5, 1e-12)

  KMK = K0*M*K0'
  inv(KMK)
  MKKMKKM = M*K0'*inv(KMK)*K0*M
  M

  JJ*MKKMKKM*JJ'

  %% Can use formula, ignoring K altogether?
  %% NO!
  F = JdaggerGtransp' * ( M ) * JdaggerGtransp;

  assert(W5, inv(F), 1e-12)

  

end
