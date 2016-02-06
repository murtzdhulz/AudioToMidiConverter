clc;
clear all;
%{
fs = 44100;
bits = 16;
recObj = audiorecorder(fs, bits, 1);

record(recObj);
disp('Start speaking')
if(strcmp(input( 'To Stop Type -1  :  '),'-1')==1)
    stop(recObj);
end
disp('End of Recording.');

%play(recObj);


myRecording = getaudiodata(recObj);
wavwrite(myRecording, fs, bits,'rec.wav');
%}
[myRecording, fs] = wavread('G_B_A_G.wav');
fft_values = fft(myRecording);
mean_value = mean(abs(fft_values));
threshold  = 1.1*mean_value;
for i=1:length(fft_values)
    if abs(fft_values(i))<threshold
        fft_values(i)=0;
    end
end
filtered_samples = ifft(fft_values);
wavwrite(filtered_samples, fs, 16,'filter.wav');
subplot(2,1,1); plot(myRecording);
subplot(2,1,2); plot(filtered_samples);
