clc;
clear all;
close all;
[y,Fs]=wavread('2 Air Choir.wav');
w=hanning(length(y));
w=w/mean(w);

left=y(:,1);
time=(1/Fs)*length(left);
t=linspace(0,time,length(left));

subplot(2,2,1);
plot(t,left);
xlabel('Time in Seconds');
ylabel('Amplitude');
title('Input Audio File in Time domain(Original)');

left=y(:,1).*w;

subplot(2,2,2);
plot(t,left);
xlabel('Time in Seconds');
ylabel('Amplitude');
title('Input Audio File in Time domain(After window is applied');

n=length(left);
p=fft(left);

nUniquePts = ceil((n+1)/2); 
p = p(1:nUniquePts); 
p = abs(p);
p=p/n;
x=p;
p=p.^2;
if rem(n, 2)
    p(2:end) = p(2:end)*2; 
  else 
    p(2:end -1) = p(2:end -1)*2; 
end


freqArray = (0:nUniquePts-1) * (Fs / n);  % create the frequency array 
subplot(2,2,4)
plot(freqArray/1000, 10*log10(p))
xlabel('Frequency (kHz)') 
ylabel('Power (dB)') 
title('Power spectrum');

subplot(2,2,3)
plot(freqArray/1000, x)
xlabel('Frequency (kHz)') 
ylabel('Magnitude') 
title('Input Audio File in Frequency domain');

[max_value, index] = max(x(:));
max_value
freqArray(index)
