function nd = fixmultirotation(dd)
%  nd = fixmultirotation(d)
% Changes d so that values are within (-pi pi)

% Kjartan Halvorsen
% 2003-12-10

  ncols = size(dd,2);
  nd=zeros(size(dd));

  for col=1:ncols
   d=dd(:,col);
   if ( max(d) > pi )
     f1 = fix( max(d) / (2*pi) );
     d = d - f1*2*pi;
   end

   if ( min(d) < - pi )
     f1 = fix( min(d) / (-2*pi) );
     d = d + f1*2*pi;
   end
  
   bigind = find( d > pi );
   for i=bigind
     d(i) = d(i) - fix(d(i)/(2*pi))*2*pi;
     if (d(i) > pi)
       d(i) = d(i) - 2*pi;
     end
   end
  
   smallind = find( d < -pi );
   for i=smallind
     d(i) = d(i) + fix(d(i)/(-2*pi))*2*pi;
     if (d(i) < -pi)
       d(i) = d(i) + 2*pi;
     end
   end
  
   nd(:,col) = d;
  end
  