function M=removeInfFFT(M)
i=1;
while(i<=length(M(:,1)))
    if M(i,3)>108 || M(i,3)<0
        M(i,:)=[];  
    end
    i=i+1;
end
end