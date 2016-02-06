load sunspot.dat
year=sunspot(:,1);
relNums=sunspot(:,2);
[pks,locs] = findpeaks(relNums);
plot(year,relNums,year(locs),pks,'rv','MarkerFaceColor','r'); grid on
xlabel('Year'); ylabel('Sunspot Number')
title('Find All Peaks'); legend('Sunspot Data','Peaks')