[y,Fs] = wavread('C_Scale.wav');

m = 2048; %no of samples in Frame
o = 512; %no of overlapping samples
n = m-o;
w = hanning(m); %hanning window
%w=w/mean(w);
left = y(:,1); %take left channel

time = (1/Fs)*length(left); %Total duration of wave file
t = linspace(0,time,length(left)); %for time axis
spectrogram(left);