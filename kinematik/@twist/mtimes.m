function t=mtimes(tw,p)
% function t=mtimes(tw,p)
% Overloading the matrix multiplication operator.
% Twist*scalar returns (4 x 4) rigid transformation, twist*point
% returns a point (3 x 1), 
% twist*twist returns a new twist.

% Kjartan Halvorsen
% 1999-05-31

if (isa(p,'twist'))
  g1=exp(tw);
  g2=exp(p);
  g=g1*g2;
  tc=expcoord(g);
  tw.coordinates=tc;
  t=twist(tw);
else
  l=length(p)
  if(l==1)
    tw.coordinates=tw.coordinates*l;
    t=gettwist(tw);
  else
    if (l==3)
      pp=[p(:);1];
    else
      pp=p;
    end
    g=exp(tw);
    tt=g*pp;
    t=tt(1:3);
  end
end
