function[odata,nmis,nmis_after]=misval_ma(idata,gapsize,order)
% function[odata,nmis,nmis_new]=misval(idata,gapsize,order)
% Program to interpolate missing data from a curve
% When the number of missing datapoints is > gapsize,
% the function is not performed. Default value for the
% acceptable gapsize is 15 frames.
% Missing data at the begin or end of the curve
% are not interpolated.
% The missing values are interpolated with a second
% order polynome (default), unless the order of the 
% polynome is given as an input variable.

%warning ('The missing values are interpolated')
mess1=['... the data has missing values at the begin of the array'];
mess2=['... the data has missing values at the end of the array'];
mess3=['... no missing values found'];
mess4=['... not enough data available for interpolating the gap'];

if nargin == 2
   order = 2;
elseif nargin == 1
   gapsize = 15;
   order = 2;
end

ndatapoints = fix(gapsize/2);

nmis=sum((idata==0));
ngood=sum((idata~=0));

[m,n]=size(idata);
odata=idata;

for i=1:n
   x=find(idata(:,i)==0);
   y=find(idata(:,i)~=0);
   
   if min(x) < min(y)
      warning(['Column ' num2str(i) ' has missing values at the begin of the array']);
   end
   
   if ~isempty(x)
      if max(x)==length(idata(:,i))
      warning(['Column ' num2str(i) ' has missing values at the end of the array']);
      end
   end
   
   %if isempty(x)
   %   warning(['Column ' num2str(i) ' has no missing values']);
   %end
   
   if length(y)>ndatapoints
      a=min(find(x>y(ndatapoints)));
   else
      a=[];
      warning(mess4);
   end
   
   while a>0
      b=min(find(y>x(a)));
      if y(b)-x(a) < gapsize
         if b+ndatapoints<=length(y) & b-ndatapoints>=1
            [V]=polyfit([y(b-ndatapoints:b+ndatapoints)],...
			idata([y(b-ndatapoints:b+ndatapoints)],i),order);
            odata(x(a):y(b)-1,i)=polyval(V,x(a):y(b)-1)';
            x=find(odata(:,i)==0);
            y=find(idata(:,i)~=0);
            a=min(find(x>y(ndatapoints)));
         else
            a=[];
         end
      else
         if isempty(b)
            a=[];
         else
            a=min(find(x>y(b)));
         end
      end
   end
   nmis_after(1,i)=sum((odata(:,i)==0));
   if ~isempty(x)
	   warning (['Column ' num2str(i) ' has ' num2str(nmis_after(1,i)) ' missing values after the interpolation procedure']);
   end
end
							      
							      


							      