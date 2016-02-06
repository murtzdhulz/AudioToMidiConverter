chirp;

x = chirp';
LENX = length(x);
sampling_rate = 1;
window = 2;
step_dist = 1;
padding = 512;
IMGY = 257;

window_length = 5;
y = STFT(x, sampling_rate, window, window_length, step_dist, padding);
AA = y';

window_length = 20;
y = STFT(x, sampling_rate, window, window_length, step_dist, padding);
BB = y';

window_length = 40;
y = STFT(x, sampling_rate, window, window_length, step_dist, padding);
CC = y';

window_length =160;
y = STFT(x, sampling_rate, window, window_length, step_dist, padding);
DD = y';


window_length = 40;
window = 1;
y = STFT(x, sampling_rate, window, window_length, step_dist, padding);
EE = y';

window = 3;
y = STFT(x, sampling_rate, window, window_length, step_dist, padding);
FF = y';

window = 4;
y = STFT(x, sampling_rate, window, window_length, step_dist, padding);
GG = y';

clf;
subplot(4,2,1)
plot(x);
xlabel('Time (seconds)');
title('Oringinal Signal');
axis([0 256 -1 1]);

freq = (1/sampling_rate)/2;

subplot(4,2,2)
imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],EE);
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title('Rectangular L = 40');
axis('xy')

subplot(4,2,3)
imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],FF);
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title('Hanning L = 40');
axis('xy')

subplot(4,2,4)
imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],GG);
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title('Blackman-Tukey L = 40');
axis('xy')

subplot(4,2,5)
imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],AA);
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title('Hamming L = 5');
axis('xy')

subplot(4,2,6)
imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],BB);
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title('Hamming L = 10');
axis('xy')

subplot(4,2,7)
imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],CC);
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title('Hamming L = 40');
axis('xy')

subplot(4,2,8)
imagesc([0:(step_dist*sampling_rate):(sampling_rate*(LENX-1))], ...
        [0:(freq/(IMGY-1)):freq],DD);
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title('Hamming L = 160');
axis('xy')







