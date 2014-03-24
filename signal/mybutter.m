function [b,a]=mybutter(n,wn)
   % Warp frequency, then find continous filter poles. Select the ones in 
   % the left half plane.
   wc=2*tan(wn*pi/2);
   k=0:2*n-1;
   args=pi/2/n*(2.*k + n -1);
   poles=wc*exp(i*args);
   ind=find(real(poles)<0);
   poles=poles(ind);

   % transform continous poles to discrete poles, using the bilinear trf
   zpoles=(2+poles)./(2-poles);
   % Compute the dc constant
   dcc=prod(1./(2-poles));

   % Compute the denominator and the nominator polynomials
   a=[1 -zpoles(1)];
   b=dcc*[1 1];
   for kk=2:length(zpoles)
      a=conv(a,[1 -zpoles(kk)]);
      b=conv(b,[1 1]);
   end
   a=real(a);
   b=real(b);
   dcgain=sum(b)/sum(a);
   b=b/dcgain;

   