[Y,p,m,S] = isp_ifchromagram('57.wav',22050);
imagesc(log(0.1*mean(Y(:))+Y))
title('Chromagram')