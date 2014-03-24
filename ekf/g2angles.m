function [theta,flag]=g2angles(jm,g)
% function [theta,flag]=g2angles(jm,g)
% Finds the solution to the inverse kinematics problem. Given the
% rigid body tranformation g, computes the solution th=[th1,..thn]
% to the problem:
%      
%     exp(w1th1)*exp(w2th2)*..*exp(wnthn)=g
%
% OBS! The method is not completely general. It is assumed that all 
% twist axes coincide at the same point. Further more, if there are 
% 6 twists in the model, the first three is for translation and the 
% last three for rotation. if the number of twists, #tw is 3<#tw<6, 
% it is assumed that the three first are for translation, the last
% one(s) for rotation. If #tw<=3 the twists are only for rotation.
%
% Input
%    jm        ->   joint model object
%    g         ->   the rigid body transformation.
%
% Output
%    theta     ->   the joint angles.
%    flag      ->   0 if ok, 1 if not
%

% Kjartan Halvorsen
% 1999-09-21

tol=1e-8;

twists=jm.twists;
ntw=length(twists);

vwq=zeros(9,ntw);
for tw=1:ntw
  vw=getCoordinates(twists{tw});
  vwq(1:6,tw)=vw;
  vwq(7:9,tw)=cross(vw(4:6),vw(1:3)); % point on the axis
end

switch ntw
  
  case 0
    % Stiff joint. Not much do do.
    theta=[];
    fl=0;
  case 1
    % Find any point not on the axis.
    [U,S,V]=svd(vwq(4:6));
    p=vwq(7:9)+U(:,3);
    p=[p;1];
    q=g*p;
    if (norm(p-q)<tol)
      theta=0;
    else
      % This is a Paden-Kahan subproblem type 1.
      [theta,fl]=padkah1(twists{1},p(1:3),q(1:3));
    end
  
  case 2
    % Find any point not on the two axes.
    w1=vwq(4:6,1);
    w2=vwq(4:6,2);
    p=vwq(7:9,2)+cross(w1,w2);
    p=[p;1];
    q=g*p;
    if (norm(p-q)<tol)
      theta=[0;0];
    else
      % This is a Paden-Kahan subproblem type 2.
      [theta,fl]=padkah2(twists{1:2},p(1:3),q(1:3));
    end
  
  case 3
    theta=zeros(3,1);
    % Find a point on the third axes that is not on the other two.
    p=intersection(twists{1:2}) + vwq(4:6,3);
    p=[p;1];
    q=g*p;
    if (norm(p-q)<tol)
      theta(1:2)=0;
    else
      % Paden-Kahan type 2är 
      [theta(1:2),fl]=padkah2(twists{1:2},p(1:3),q(1:3));
    end
    
    % Solve for the remaining angle
    % A point not on the third axis
    [U,S,V]=svd(vwq(4:6,3));
    p=vwq(7:9,3)+U(:,3);
    p=[p;1];
    g1g2=trf(twists{1},theta(1))*trf(twists{2},theta(2));
    q=inv(g1g2)*g*p;
    if (norm(p-q)<tol)
      theta(3)=0;
    else
      % Paden-Kahan type 1
      [theta(3),fl]=padkah1(twists{3},p(1:3),q(1:3));
    end
    
  case 4
    theta=zeros(4,1);
    % Find the point of intersection of the last twists.
    p=intersection(twists{2:3});
    p=[p;1];
    q=g*p;
    dp=q(1:3)-p(1:3);
    theta(1)=(dp'*vwq(1:3,1));
    g1=trf(twists{1},theta(1));
    gg=inv(g1)*g;
    % Find the next three angles as above
    % Find a point on the third axes that is not on the other
    % two. p from above
    p=p(1:3) + vwq(4:6,4);
    p=[p;1];
    q=gg*p;
    if (norm(p-q)<tol)
      theta(2:3)=0;
    else
      % Paden-Kahan type 2 
      [theta(2:3),fl]=padkah2(twists{2:3},p(1:3),q(1:3));
    end
    
    % Solve for the remaining angle
    % A point not on the third axis
    [U,S,V]=svd(vwq(4:6,4));
    p=vwq(7:9,4)+U(:,3);
    p=[p;1];
    g2g3=trf(twists{2},theta(2))*trf(twists{3},theta(3));
    q=inv(g1g2)*gg*p;
    if (norm(p-q)<tol)
      theta(4)=0;
    else    
      % Paden-Kahan type 1
      [theta(4),fl]=padkah1(twists{4},p(1:3),q(1:3));
    end
  
  case 5
    theta=zeros(5,1);
    % Find the point of intersection of the last twists.
    p=intersection(twists{3:4});
    p=[p;1];
    q=g*p;
    dp=q(1:3)-p(1:3);
    theta(1:2)=(dp'*vwq(1:3,1:2))';
    g1g2=trf(twists{1},theta(1))*trf(twists{2},theta(2));
    gg=inv(g1g2)*g;
    % Find the next three angles as above
    % Find a point on the third axes that is not on the other two.
    p=p(1:3) + vwq(4:6,5);
    p=[p;1];
    q=gg*p;
    if (norm(p-q)<tol)
      theta(3:4)=0;
    else
      % Paden-Kahan type 2 
      [theta(3:4),fl]=padkah2(twists{3:4},p(1:3),q(1:3));
    end
    
    % Solve for the remaining angle
    % A point not on the third axis
    [U,S,V]=svd(vwq(4:6,5));
    p=vwq(7:9,5)+U(:,3);
    p=[p;1];
    g3g4=trf(twists{3},theta(3))*trf(twists{4},theta(4));
    q=inv(g3g4)*gg*p;
    if (norm(p-q)<tol)
      theta(5)=0;
    else
      % Paden-Kahan type 1
      [theta(5),fl]=padkah1(twists{5},p(1:3),q(1:3));
    end
    
  case 6
    theta=zeros(6,1);
    % Find the point of intersection of the last twists.
    p=intersection(twists{4:5});
    p=[p;1];
    q=g*p;
    dp=q(1:3)-p(1:3);
    vwq(1:3,1:3);
    theta(1:3)=(dp'*vwq(1:3,1:3))';
    g1g2g3 = ...
	trf(twists{1},theta(1)) * trf(twists{2},theta(2)) * ...
	trf(twists{3},theta(3));
    gg=inv(g1g2g3)*g;
    % Find the next three angles as above
    % Find a point on the third axes that is not on the other two.
    p=p(1:3) + vwq(4:6,6);
    p=[p;1];
    q=gg*p;
    if (norm(p-q)<tol)
      theta(4:5)=0;
    else
      % Paden-Kahan type 2 
      [theta(4:5),fl]=padkah2(twists{4:5},p(1:3),q(1:3));
    end
    
    % Solve for the remaining angle
    % A point not on the third axis
    [U,S,V]=svd(vwq(4:6,6));
    p=vwq(7:9,6)+U(:,3);
    p=[p;1];
    g4g5=trf(twists{4},theta(4))*trf(twists{5},theta(5));
    q=inv(g4g5)*gg*p;

    if (norm(p-q)<tol)
      theta(6)=0;
    else
      % Paden-Kahan type 1
      [theta(6),fl]=padkah1(twists{6},p(1:3),q(1:3));
    end
    
end

flag=fl;

  
  
