function [M,notearray]=detectPoly(left,m,n,Fs,w,k,thrfact,amtpoly)
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
    [pks,locs]=findpeaks(p,'SORTSTR','descend');
    notearray(i,1:amtpoly)=freqArray(locs(1:amtpoly));
    maxvalarr(i,1:amtpoly)=pks(1:amtpoly);
    i = i+1;
    s1 = s1+n-1;
end

thr=thrfact*max(maxvalarr);

for i=2:k-1
    for j=1:amtpoly
        if(maxvalarr(i,j)<thr(j))
            maxvalarr(i,j)=maxvalarr(i-1,j);
            notearray(i,j)=notearray(i-1,j);
        end
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

x=1;
for j=1:amtpoly
    st = 0; %Start Time
    i=1;
    while i < k
        M(x,3)=round(12*log2(notearray(i,j)/440)+69);
        M(x,4)=(1+mapminmax(maxvalarr(i,j)))*60;
        nt=0; %note time
        while i<k && (notearray(i,j)==notearray(i+1,j)) 
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
end
M(:,1) = 1;         % all in track 1
M(:,2) = 1;         % all in channel 1
%M(:,4) = 100;  % lets have volume ramp up 80->120
end