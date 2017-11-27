function  [IMall, indAcc, Mall, gcnames,accAll] = interaction_moments_base_two_arms(gmbase, gmleft, gmright, ...
                                        statesbase, statesleft, statesright,...
                                        accBase, accLeft, accRight)
     % Combines the interaction moments of the base and two arms
     % Input
     %    gmbase, gmleft, gmright   ->  Models for the base (hip and
     %                                  trunk), the left chain and right
     %                                  chain.
     %    statesbase, statesleft, statesright  -> (2*ndofs x nfrs) state
     %                                  vectors with joint angles and
     %                                  angular velocities.
     % Output
     %    IMall       <-  The interaction moments (ndofs x nfrs)
     %    IndAcc      <-  The accelerations induced by the interaction moments (ndofs x nfrs)
     %    Mall        <-  The generalized manipulator inertia (ndofs x ndofs x nfrs)
     %    gcnames     <-  Cell array with names of the joint angles /
     %                    degress of freedom


     nfrs = size(statesleft, 2);
     nBaseStates = size(statesbase, 1)/2;
     nLeftArmStates = size(statesleft, 1)/2 - nBaseStates;
     nRightArmStates = size(statesright, 1)/2 - nBaseStates;
     nBS = nBaseStates;
     nLS = nLeftArmStates;
     nRS = nRightArmStates;
    
     Mbase = generalized_manipulator_inertia(gmbase, statesbase(1:nBS, :));
     Mleft = generalized_manipulator_inertia(gmleft, statesleft(1:nBS+nLS, :));
     Mright = generalized_manipulator_inertia(gmright, statesright(1:nBS+nRS, :));

     % Debug. OK!
     %[twsBase, g0Base, MbBase, gcnamesBase] = flatten_km(gmbase);
     %[twsLeft, g0Left, MbLeft, gcnamesLeft] = flatten_km(gmleft);
     %[twsRight, g0Right, MbRight, gcnamesRight] = flatten_km(gmright);
     %Mbase_flat = generalized_manipulator_inertia_flattened(MbBase, twsBase, g0Base, statesbase(1:nBS, 60));
     %Mleft_flat = generalized_manipulator_inertia_flattened(MbLeft, twsLeft, g0Left, statesleft(1:nLS+nBS, 60));
     %Mright_flat = generalized_manipulator_inertia_flattened(MbRight, twsRight, g0Right, statesright(1:nRS+nBS, 60));
     
     Mbase = Mbase(1:nBS, 1:nBS, :); % Kolla om inte har nBS rader från början
     Mleft = Mleft(1:nBS+nLS, 1:nBS+nLS, :);
     Mright = Mright(1:nBS+nRS, 1:nBS+nRS, :);
     
     
     % Combine
     ML_B = Mleft(1:nBS, 1:nBS,:);
     ML_L = Mleft(nBS+1:end,nBS+1:end,:);
     ML_BL = Mleft(1:nBS, nBS+1:end,:);
     MR_B = Mright(1:nBS, 1:nBS,:);
     MR_R = Mright(nBS+1:end,nBS+1:end,:);
     MR_BR = Mright(1:nBS, nBS+1:end,:);
 
     Mall = zeros(nBS+nLS+nRS, nBS+nLS+nRS, nfrs);
     Mall(1:nBS, 1:nBS, :) = Mbase+ML_B+MR_B;
    
     Mall(1:nBS, nBS+1:nBS+nLS, :) = ML_BL;
     Mall(nBS+1:nBS+nLS, 1:nBS,:) = permute(ML_BL, [2,1,3]);
     Mall(nBS+1:nBS+nLS, nBS+1:nBS+nLS,:) = ML_L;
 
     Mall(1:nBS, nBS+nLS+1:end, :) = MR_BR;
     Mall(nBS+1+nLS:end, 1:nBS,:) = permute(MR_BR, [2,1,3]);
     Mall(nBS+1+nLS:end, nBS+1+nLS:end,:) = MR_R;
     
 
     accAll = cat(1, accBase, accLeft(nBS+1:end,:), accRight(nBS+1:end, :));
 
     % Calculate the interaction due to acceleration at other angles (degrees of
     % freedom). Remove the diagonal, then multiply Mall with acc vector
     ntot = nBS+nLS+nRS;
     IM_acc = zeros(ntot, nfrs);
     for i=1:nfrs
         Mi = Mall(:,:,i); 
         Mi_nodiagonal = Mi - diag(diag(Mi));
         % IM_acc(:,i) = - Mi_nodiagonal*cat(1, accBase(:,i), accLeft(nBS+1:end,i), accRight(nBS+1:end,i));
         IM_acc(:,i) = - Mi_nodiagonal*accAll(:,i);
     end
     

     %% Calculate the interaction moments. OBS these are due to velocity only. Must be combined with IM_acc
     % im is the interaction moment for each degree of freedom, C
     % is the Coriolis matrix, i.e. 
     %  im = -C * \dot{q} 
     disp('Computing interaction moments for left arm')
     [IMleft, Cleft, dofnamesLeft] = interaction_moments(gmleft, statesleft); %, dofs2study);
     disp('Computing interaction moments for right arm')
     [IMright, Cright, dofnamesRight] = interaction_moments(gmright, statesright);%, dofs2study);
     disp('Computing interaction moments for hip and trunk')
     [IMbase, Cbase, dofnamesBase] = interaction_moments(gmbase, statesbase);%, dofs2study);

     %% Combine interaction terms 
     
     
     IMbaseLR = IMbase + IMleft(1:nBaseStates, :) ...
         + IMright(1:nBaseStates, :);

     IM_vel = cat(1, IMbaseLR, IMleft(nBaseStates+1:end,:), ...
         IMright(nBaseStates+1:end, :));
     IMall =  IM_vel + IM_acc;
     
 
     % When calculating the Coriolis matrix, we assume the vector
     % of joint velocities to be [\dot{q}_B, \dot{q}_L, \dot{q}_R]
     %CbaseLR = zeros(nBaseStates, ...
     %                nBaseStates+nLeftArmStates+nRightArmStates, ... 
     %                nfrs); 
     %CbaseLR(:,1:nBaseStates, :) = Cbase ...
     %    + Cleft(1:nBaseStates, 1:nBaseStates,:) ...
     %    + Cright(1:nBaseStates, 1:nBaseStates,:);
     %CbaseLR(:, nBaseStates+1:nBaseStates+nLeftArmStates,:) = ...
     %    Cleft(1:nBaseStates, nBaseStates+1:end, :);
     %CbaseLR(:, nBaseStates+nLeftArmStates+1:end, :) = ...
     %        Cright(1:nBaseStates, nBaseStates+1:end, :);
             

     %% Now the induced accelerations. Since we have taken care of the off-diagonal 
     % inertial terms these should not be included
     indAcc = zeros(ntot, nfrs);
     for i = 1:nfrs
         indAcc(:,i) = IMall(:, i)./diag(Mall(:,:,i));
     end
     
     %% Finally, the names
     gcnames = cat(1, dofnamesBase, dofnamesLeft(nBS+1:end), dofnamesRight(nBS+1:end));
 
