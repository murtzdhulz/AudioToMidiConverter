clc;
clear all;
close all;
[y,Fs] = wavread('Piano.wav');
%{
N = length(y);
totaldur = length(y)/Fs;
t = linspace(0,totaldur,N);
n = y(t>0.9 & t<1.9); 
n = vertcat(n,n,n,n,n,n,n,n,n,n,n,n,n);
y=noiseRemove(y,n);
%}
m = 2048; %no of samples in Frame
o = 512; %no of overlapping samples
n = m-o;
w = hanning(m); %hanning window
%w=w/mean(w);
left = y(:,1); %take left channel

time = (1/Fs)*length(left); %Total duration of wave file
t = linspace(0,time,length(left)); %for time axis

subplot(2,1,1);
plot(t,left);
xlabel('Time in Seconds');
ylabel('Amplitude');
title('Input Audio File in Time domain');

thr=0.01;
left=clip(left,thr);


subplot(2,1,2);
plot(t,left);
xlabel('Time in Seconds');
ylabel('Amplitude');
title('Input Audio File in Time domain');

[left, k]=padAudio(left,m,n);

[M,notearray]=detectPoly(left,m,n,Fs,w,k,0.2,2);

M=removeInfFFT(M);

M=equalizeNotes(M,n,Fs,2);

M=removeInfFFT(M);
%%
%{
[pks,locs] = findpeaks(abs(left));

figure(100);
plot(locs,pks);
%}
%%
%polySynthAdd(M,Fs,0.3,0.6,'sin','sin','square',0.7,1,0.1,0.8,0.1,0.1);
%polySynthFM(M,Fs,1,5,'sin','sin',0.7,1,0.1,0.8,0.1,0.1);
%polyClarinet(M,Fs,0.1,1,0.3,0.8,0.4,0.2);
%%
midi_new = matrix2midi(M);
writemidi(midi_new, 'testout.mid');
%soundsc(midi2audio('testout.mid',Fs,'fm'),Fs);
%%
ispmidi = isp_midiread('testout.mid');
ispmidi.instruments=zeros(1,4);
ispmidi.controller=zeros(1,5);
[wav, fs]=isp_midisynth(ispmidi, 'soundfont', 'stavi_violin.sf2');
soundsc(wav,fs);
%%
figure(2);
plot(notearray);
xlabel('Frame number'); 
ylabel('Frequency');
%% convert to 'Notes' matrix:
Notes = midiInfo(midi_new,0);

%% compute piano-roll:
[PR,t,nn] = piano_roll(Notes,1);

%% display piano-roll:
figure;
imagesc(t,nn,PR);
axis xy;
xlabel('time (sec)');
ylabel('note number');

figure;
nmat = midi2nmat('testout.mid');
pianoroll(nmat,'Test','sec','vel');
%playsound(nmat,'shepard',Fs);
%playmidi(nmat);
%%
%figure;
%[u,v,w]=myspecgram(left,m,Fs,[],o);
%spectrogram(left);