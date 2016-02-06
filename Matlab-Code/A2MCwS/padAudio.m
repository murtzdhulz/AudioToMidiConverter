function [left, k]=padAudio(left,m,n)
al = length(left)-m;
pad1 = mod(al,n);
k = ceil(al/n) + 1;

left = padarray(left,[pad1+m 0],0,'post');
end