clc;
clear all;
close all;
[y,Fs]=wavread('audiocheck.net_sin_440Hz_-3dBFS_1s.wav');
%plot(y);
R=10;
w=hanning(length(y));
w=w/mean(w);
y=y(:,2);
yw=w.*y;
plot(yw);
wavefft=abs(fft(yw,1000));
figure(2);
plot(wavefft);
w=hamming(length(y));
xlabel('Frequency in Hz');
ylabel('Magnitude');
title('The Wave FFT');

spectrogram(y,1024,128,256,Fs); 
figure(4);
nfft = 2^(nextpow2(length(yw)));   % Find next power of 2
fftx = abs(fft(yw,nfft));
plot(fftx);

figure(5);
Hs1 = spectrum.mtm(4,'adapt');
psd(Hs1,yw,'Fs',Fs,'NFFT',1024)
