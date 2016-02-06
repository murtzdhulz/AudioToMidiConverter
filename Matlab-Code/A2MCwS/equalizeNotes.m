function M=equalizeNotes(M,n,Fs,x)
i=1;
while(i<length(M(:,1)))
    if M(i,6)-M(i,5)<=x*(n/Fs)
        M(i+1,5)=M(i,5);
        M(i,:)=[];
    end
    i=i+1;
end
for i=2:length(M(:,3))
    if (M(i-1,3)-M(i,3))>10 
        M(i,3)=M(i,3)+12;
    end
    if (M(i-1,3)-M(i,3))<-10
        M(i,3)=M(i,3)-12;
    end
end
%{
for i=1:length(M(:,3))-1
    if (M(i,3)-M(i+1,3))>10 
        M(i,3)=M(i,3)-12;
    end
    if (M(i,3)-M(i+1,3))<-10
        M(i,3)=M(i,3)+12;
    end
end
%}

M=sortrows(M,5);
i=1;
while(i<length(M(:,3)))
    if M(i,3)==M(i+1,3)
        M(i,6)=M(i+1,6);
        M(i+1,:)=[];
    end
    i=i+1;
end

end