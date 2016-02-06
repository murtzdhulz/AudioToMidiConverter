%ISP_CHROMA  Compute chroma vectors
%
% SYNTAX
%   [vc, options]= isp_chroma(y, options ...)
%
% DESCRIPTION
%   Computes the chroma vectors as given by Goto [Masataka Goto: A
%   Chorus-Selecting Detecting Method for Musical Audio Signals]. The
%   idea is that each vector corresponds to a note (a,c,f etc).
%
% INPUT
%   y:
%     Wave signal.
%   options ...:
%     Structs or field/value pairs specifying any of the following:
%     octL:
%       Lowest octave (default: 3).
%     octH:
%       Highest octave (default: 8).
%     sampf:
%       Sampling frequency (default: 44100).
%     windowsize:
%       Step size (default: 80ms).
%     fftwindow:
%       FFT window size. The frame is modified to match 2^x (default: 250ms).
%
% OUTPUT
%   vc:
%     12 dimensional chroma vector
%   options:
%     Input arguments supplemented by default values.
%
% SEE ALSO
%   isp_ifchromagram.
%
% HISTORY
%   21-04-2006: Created by Mads Emil Solgård, Filip, Tue Lehn-Schiøler.
%   2008: Modified to toolbox format by Jesper Højvang Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [vc, par]= isp_chroma(y, varargin)

par = isp_interpretarguments(struct( ...
    'octL', 3, ...
    'octH', 8, ...
    'sampf', 44100, ...
    'windowsize', 80, ...
    'fftwindow', 250, ...
    'timestamp', datestr(now)), varargin{:})


fs=par.sampf;
fw = par.fftwindow;
%Denne del udregner Short-time Fourier Transform (STFT)

shift = fs*par.windowsize/1000;      % shifts that are default 80ms long
wav_length = length(y);
length_segment = 2^floor(log2(fs*fw/1000));

num_segment = round(wav_length/shift-length_segment/shift)    ;

win = hanning(length_segment)';

for i= 1:num_segment
    data(i,:)=y( (i-1)*shift+1 : (i-1)*shift+length_segment);
    data(i,:) = data(i,:).*win;
    spec(i,:) = abs(fft(data(i,:))).^2;
end



%Denne del "konstruerer" De Band pass filtre der skal bruges

hertz_bin = (fs/2)/(length_segment/2);

BPF = zeros(12, length_segment/2);

for c = 1:12
    for h = par.octL:par.octH
        Fch = 1200*h+100*(c-1);

        fcent = 1200*h+100*(c-1-1);
        f = 2^(fcent/1200)*440*2^(3/12-5);
        bin_number = ceil(f/hertz_bin);

        fcentH = 1200*h+100*(c);
        fH = 2^(fcentH/1200)*440*2^(3/12-5);
        bin_numberH = floor(fH/hertz_bin);

        while bin_number < bin_numberH

            BPF(c, bin_number) = 1/2*(1-cos(2*pi*(fcent-(Fch-100))/200));

            bin_number = bin_number +1;
            f = f + hertz_bin;
            fcent = 1200*log2(f/(440*2^(3/12-5)));
        end;
    end;
end;

%Udregn chroma vektorerne ved at multiplicere filteret med spektret.
vc = BPF * spec(:, 1:length_segment/2)';

% figure;
% imagesc(vc)

%save speak4min vc
