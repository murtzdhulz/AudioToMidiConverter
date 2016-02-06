function left=clip(left,thr)
for i=1:length(left)
    if abs(left(i))<thr
        left(i)=0;
    end
end
end