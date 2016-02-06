clear('all');
mu;
mu1;
mu2;

mu = fliplr(mu');
mu1 = fliplr(mu1');
mu2 = fliplr(mu2');

mu2 = log(mu2);
MUX = mu2-mean(mu2);

LENX = length(MUX);
sampling_rate = 1;
window = 2;
step_dist = 1;
padding = 512;
IMGY = 257;

window_length = 100;
y = STFT(MUX, sampling_rate, window, window_length, step_dist, padding);
AA = y';

clf;
subplot(4,1,1)
plot(mu);
xlabel('Time (0.5 months)');
title('MU Original');
axis([1 256 0 90]);

subplot(4,1,2)
plot(MUX);
xlabel('Time (0.5 months)');
title('MU Not Split Adjusted, Dollar Cost Factored, Log, Zero Mean');
axis([1 256 -1.5 1.5]);

freq = (1/sampling_rate)/2;

subplot(4,1,3)
imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],AA);
xlabel('Time (0.5 months)');
ylabel('Normalized Frequency');
title('STFT of Above, Hamming Window Length = 100');
axis('xy')

subplot(4,1,4)
BB = log(AA);
BB = BB - min(min(BB));
BB = (BB/(max(max(BB))))*255;
BB(13,:) = 265*ones(1,256);
BB(5,:) = 265*ones(1,256);

imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],BB);
xlabel('Time (0.5 months)');
ylabel('Normalized Frequency');
title('Log and Zoom of Above with Markers Added');
axis([0 256 0 0.1]) 
axis('xy')
