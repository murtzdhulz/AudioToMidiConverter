%ISP_SHOWSPECTROGRAM  Plots a spectrogram
%
% SYNTAX
%     [freq, time] = isp_showspectrogram(X, options ...)
%
% DESCRIPTION
%   Plots a spectrogram with correct axes.
%
% INPUT
%   X:
%     The spectrogram to be plotted, e.g. obtained with ...
%   options ...:
%     Structs and 'field'/value pairs:
%     samplerate:
%       Sample frequency of X.
%     fftlength:
%       FFT length.
%     hopsize:
%       The window hop size. Default: half the FFT length.
%
% OUTPUT
%   freq:
%     Tick values on frequency axis.
%   time:
%     Tick values on time axis.
%
% HISTORY
%   Created by Jesper H. Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function isp_showspectrogram(X, varargin)
    options.samplerate=1;
    options.hopsize=0;
    options.fftlength=0;
    
    options = isp_interpretarguments(options, varargin{:});
    if options.fftlength==0
        options.fftlength = size(X, 1);
    end

    if options.hopsize==0
        options.hopsize = options.fftlength/2;
    end

    freq=(0:size(X, 1)-1)*options.samplerate/options.fftlength;
    time=(0:size(X, 2)-1)*options.hopsize/options.samplerate;

    % Don't include stuff above half the sample rate
    uniqFreq = freq<=options.samplerate/2;

    surf(time, freq(uniqFreq), X(uniqFreq,:), 'EdgeColor', 'none')

    view(0,90)
    axis tight
    axis xy

    set(gca, 'YScale', 'log')
    xlabel('Time')
    ylabel('Frequency')
    title('Spectrogram')
end