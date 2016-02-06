clc;
clear all;
close all;
[y,Fs] = wavread('GuitarP.wav');

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

subplot(2,1,2);
plot([1:m],w);
xlabel('No. of samples');
ylabel('Weighting');
title(['Hanning Window of ' num2str(m) ' Samples']);

al = length(left)-m;
pad1 = mod(al,n);
k = ceil(al/n) + 1;

left = padarray(left,[pad1+m 0],0,'post');

nUniquePts = ceil((m+1)/2);
freqArray = (0:nUniquePts-1)*(Fs/m);

% build matrix of windowed data slices
ftime = n/Fs;

s1 = 1;
i = 1;
notearray = zeros(1,k);

while i <= k
    frame = left(s1:(s1+m-1));
    wframe = frame.*w;
    fl = length(frame);
    p = fft(wframe);
    p = p(1:nUniquePts);
    p = abs(p);
    p = p/fl;
    p = p.^2;
    if rem(fl, 2)
        p(2:end) = p(2:end)*2;
    else
        p(2:end-1) = p(2:end-1)*2;
    end
    [max_value, index] = max(p(:));
    [c,index1] = max(p(p~=max(p(:))));
    if index1>index
        index1=index1+1;
    end
    notearray(i) = freqArray(index);
    notearray1(i)= freqArray(index1);
    i = i+1;
    s1 = s1+n-1;
end
%%
%{
M = zeros(k,6);
M(:,1) = 1;         % all in track 1
M(:,2) = 1;         % all in channel 1
M(:,3) = 12*log2(notearray/440)+69;      % note numbers: one ocatave starting at middle C (60)
M(:,4) = 100;  % lets have volume ramp up 80->120
M(:,5) = (0:(n/Fs):time)';  % note on:  notes start every .5 seconds
M(:,6) = M(:,5) + (n/Fs);   % note off: each note has duration .5 seconds
%}

st = 0; %Start Time
x=1;
i=1;
while i < k
    M(x,3)=round(12*log2(notearray(i)/440)+69);
    nt=0; %note time
    while i<k && (notearray(i)==notearray(i+1)) 
        i=i+1;
        nt = nt+ftime;
    end
    nt = nt+ftime;
    M(x,5)=st;
    M(x,6)=st+nt;
    st=st+nt;
    x=x+1;
    i=i+1;
end
st = 0;
while i < 2*k
    M(x,3)=round(12*log2(notearray1(i-k)/440)+69);
    nt=0; %note time
    while i<k && (notearray1(i-k)==notearray1(i-k+1)) 
        i=i+1;
        nt = nt+ftime;
    end
    nt = nt+ftime;
    M(x,5)=st;
    M(x,6)=st+nt;
    st=st+nt;
    x=x+1;
    i=i+1;
end
M(:,1) = 1;         % all in track 1
M(:,2) = 1;         % all in channel 1
M(:,4) = 100;  % lets have volume ramp up 80->120

for i=2:length(M(:,3))-1
    if (M(i-1,3)-M(i,3))>10 
        M(i,3)=M(i,3)+12;
    end
    if (M(i-1,3)-M(i,3))<-10
        M(i,3)=M(i,3)-12;
    end
end
midi_new = matrix2midi(M);
writemidi(midi_new, 'testout.mid');
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
playsound(nmat,'fm',Fs);
%%
%figure;
%[u,v,w]=myspecgram(left,m,Fs,[],o);
%spectrogram(left);