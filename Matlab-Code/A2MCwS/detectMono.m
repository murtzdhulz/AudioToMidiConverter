function [M,notearray]=detectMono(left,m,n,Fs,w,k,thrfact)
nUniquePts = ceil((m+1)/2);
freqArray = (0:nUniquePts-1)*(Fs/m);

% build matrix of windowed data slices
ftime = n/Fs;

s1 = 1;
i = 1;
notearray = zeros(1,k);
maxvalarr = zeros(1,k);
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
    maxvalarr(i) = max(frame);
    notearray(i) = freqArray(index);
    i = i+1;
    s1 = s1+n-1;
end

thr=thrfact*max(maxvalarr);

for i=2:k-1
    if(maxvalarr(i)<thr)
        maxvalarr(i)=maxvalarr(i-1);
        notearray(i)=notearray(i-1);
    end
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
    M(x,4)=(1+mapminmax(maxvalarr(i)))*60;
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
M(:,1) = 1;         % all in track 1
M(:,2) = 1;         % all in channel 1
%M(:,4) = 100;  % lets have volume ramp up 80->120
end