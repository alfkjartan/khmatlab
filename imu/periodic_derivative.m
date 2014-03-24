function dx = periodic_derivative(x, dt, plotit)
  %% Derivative using fft

  if (nargin == 0)
    demo_me();
    return
  end

  %% Calculate fourier transform
  %% fftshift reorders the resulting vector to represent negative
  %% frequencies in the first half of the vector
  [N, dim] = size(x);

  f = ((1:N) - ceil(N/2))' / N / dt;
  %%keyboard
  %%f = -(nfr/2):1:(nfr/2); 
  %% Time derivative in the frequency domain
  wi=2*pi*f*i;

  dx = zeros(size(x));
  for dd = 1:dim
    b=fftshift(fft(x(:,dd)));
    c=b.*wi;

 % Inverse FFT
 % Although the result should be real ideally, there will be
 % small imaginery parts due to numeric inaccuracies, thus plotting
 % only the real part
    d=ifft(ifftshift(c));
    dx(:,dd) = real(d);
  end

 if (plotit)
 subplot(4,1,1), plot(x)
 subplot(4,1,2), plot(f,abs(b));
 subplot(4,1,3), plot(f,abs(c));
 subplot(4,1,4), plot(real(d));
end

function demo_me()
  t=0:0.01:0.99;
  a=sin(2*pi*10*t);

  da = periodic_derivative(a', 0.01, 1);
