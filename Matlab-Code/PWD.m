function y = pwd(x,fftl)

%
%  y = pwd(x,fftl)
%
%  pwd produces a Pseudo-Wigner Distribution of "x".
%  The output is stored in "y".
%
%  x = input signal to be analyzed
%  fftl = Length of FFT to be used
%
%  Developed by Srikrishna Bhashyam
%               Rice University
%               May, 1999
%               skrishna@rice.edu
%
%  Coded using MATLAB 5.X.X.
%
%	REVISION HISTORY
%
%	VERSION 1.0.0		MAY 5, 1999	Srikrishna Bhashyam
%

L = length(x);

wlength = L;
%  wlength = truncation length of wigner dist to get pseudo-wigner dist

zp = fftl - 2*L + 1;

y = zeros(fftl,L);
xx = zeros(1,2*wlength + 1);

xinter = interp(x,2);
xinter = xinter(1:2*L-1);

for k = 1:L
  if ((-wlength/2 + k) < 1)
    xx1 = [zeros(1,wlength - 2*k) xinter(1:2*(wlength/2 + k)-1)];
  else
    xx1 = [xinter(2*(-wlength/2+k)-1:2*L-1) zeros(1,2*k-2-wlength)];
  end
  xx2 = conj(fliplr(xx1));
  xx = xx1.*xx2;
  xxzp = [xx(L:2*L-1) zeros(1,zp) xx(1:L-1)];
%  yktemp = (fft(xxzp).*exp(-sqrt(-1)*2*pi*[0:fftl-1]*wlength/fftl))/L;
  yktemp = fft(xxzp)/sqrt(2*L-1);
  y(:,k) = conj(yktemp(1:fftl)');
end

imagesc([1:L],[0:fftl-1]/fftl,real(y));
axis('xy');



