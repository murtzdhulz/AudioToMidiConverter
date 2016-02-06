clc; close all; clear all;
% didn't have wav file, but simply replace this with the following
 [audio,fs] = wavread('57.wav');
 audio= audio(:,1);
%audio = rand(1,10000);
fs = 44100; % temp sampling frequency, will depend on audio input
NFFT = 1024; % feel free to change FFT size
hamWin = hamming(NFFT); % window your audio signal to avoid fft edge effects

% get spectral content
S = spectrogram(audio,hamWin,NFFT/2,NFFT,fs);

% Start at center lowest piano note
A0 = 27.5;
% all 88 keys
keys = 0:87;
center = A0*2.^((keys)/12); % set filter center frequencies
left = A0*2.^((keys-1)/12); % define left frequency
left = (left+center)/2.0;
right = A0*2.^((keys+1)/12); % define right frequency
right = (right+center)/2;

% Construct a filter bank
filter = zeros(numel(center),NFFT/2+1); % place holder
freqs = linspace(0,fs/2,NFFT/2+1); % array of frequencies in spectrogram
for i = 1:numel(center)
    xTemp = [0,left(i),center(i),right(i),fs/2]; % create points for filter bounds
    yTemp = [0,0,1,0,0]; % set magnitudes at each filter point
    filter(i,:) = interp1(xTemp,yTemp,freqs); % use interpolation to get values for   frequencies
end

% multiply filter by spectrogram to get chroma values.
chroma = filter*abs(S);

%Put into 12 bin chroma
chroma12 = zeros(12,size(chroma,2));
for i = 1:size(chroma,1)
    bin = mod(i,12)+1; % get modded index
    chroma12(bin,:) = chroma12(bin,:) + chroma(i,:); % add octaves together
end

plot(chroma12);