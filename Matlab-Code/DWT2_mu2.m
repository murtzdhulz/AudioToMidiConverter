clear('all');
mu;
mu1;
mu2;

mu = fliplr(mu');
mu1 = fliplr(mu1');
mu2 = fliplr(mu2');

mu2 = log(mu2);
MUX = mu2-mean(mu2);
sampling_rate = 1;

clf;
subplot(3,1,1)
scale_level =  4;
y = DWT2(MUX, sampling_rate, scale_level);
xlabel('Time (0.5 months)');
ylabel('Normalized Frequency');
title('Scale Level = 4');
axis([0 255 0 0.5]);

subplot(3,1,2)
scale_level =  8;
y = DWT2(MUX, sampling_rate, scale_level);
xlabel('Time (0.5 months)');
ylabel('Normalized Frequency');
title('Scale Level = 8');
axis([0 255 0 0.5]);

subplot(3,1,3)
scale_level =  8;
y = DWT2(MUX, sampling_rate, scale_level);
xlabel('Time (0.5 months)');
ylabel('Normalized Frequency');
title('Scale Level = 8, Zoom of Above');
axis([0 255 0 0.1]);
