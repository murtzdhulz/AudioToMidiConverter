%ISP_WAVMODIFY  Modifies properties of a WAV struct.
%
% SYNTAX
%   [wavout, optsout] = isp_wavmodify(wavin, options, 'field', value, ...)
%   
% DESCRIPTION
%   Add silence to a signal or modify its the bitrate, SNR or bandwidth
%   of a signal.
%
% INPUT
%   wavin:
%     Structure describing a wav song.
%   options, field/value pairs:
%     The following parameters can be set as field names in the
%     'options' struct or be specified as field/value pairs:
%     removesilence:
%       Remove silence at the beginning and end of the signal.
%     samplerate:
%       Input sample rate. Required for bitrate and bandwidth
%       modifications.
%     snr:
%       Add white, Gaussian noise with the SNR specified in dB. Default: Inf.
%     bitrate:
%       Compress and decompress with the lame mp3 encoder at the bitrate
%       given in kbps to introduce compression artifacts. Inf denotes no
%       compression. Default: Inf
%     mono:
%       Average all channels into a mono signal. Default: false.
%     bandwidth:
%       Remove all frequency information above 'bandwidth'. Requires
%       'samplerate' to be specified as well. Default: Inf.
%     maxlength:
%       The maximum length in seconds of the song. If the song is
%       longer, only the middle part is retained. Default: Inf.
%     addsilence:
%       Specifies how much silence to add to the signal. It is specified
%       as the fraction of the length of the original signal. Noise is
%       added on the order of the least significant bit. Default: 0.
%
% OUTPUT
%   wavout:
%     wav song. The modifications are applied in the order mentioned under
%     'options, field/value pairs'. Thus, specifying a maxlength of 100 s
%     and 'addsilence' 0.1, the resulting signal will be of length 110 s.
%   optsout:
%     Structure specifying the actual modification parameters.
%
% SEE ALSO
%   isp_readwav, isp_writewav, isp_wavdemo.
%   
% HISTORY
%   Created by Jesper H. Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [wav, options]=isp_wavmodify(wav, varargin)
    options = struct('removesilence', true, ...
                     'samplerate', nan, ...
                     'bitrate', inf, ...
                     'snr', inf, ...
                     'bandwidth', inf, ...
                     'mono', false, ...
                     'maxlength', inf, ...
                     'addsilence', 0);
    options = isp_interpretarguments(options, varargin{:});

    if isnan(options.samplerate) && ...
            (~isinf(options.bitrate) || ~isinf(options.bandwidth))
        error('Samplerate not specified.')
    end

    if options.removesilence

        % Remove silence at the beginning and in the end
        % A for-loop seems to be quicker than variants of 
        %   find( all(abs(wav)>=eps+1/32768, 2) )
        % Even when using the 'first' or 'last' options, the inequality
        % still has to be evaluated for the entire 'wav'.
        for wBeg=1:length(wav)
            if any(abs(wav(wBeg,:)) >= eps+1/32768 )
                break
            end
        end
        for wEnd=length(wav):-1:1
            if any(abs(wav(wEnd,:)) >= eps+1/32768 )
                break
            end
        end

        wav=wav(wBeg:wEnd,:);
    end
          
    % Change SNR
    if ~isinf(options.bitrate)
        wavStd = std(wav(:));
        noiseStd = wavStd / 10^(options.snr/20);
        wav = wav + noiseStd*randn(size(wav));
    end

    % Change bitrate
    if ~isinf(options.bitrate)
        nCh = size(wav, 2);
        if  nCh ~= 2
            warning(['Input is not stereo. Bitrates might not be ' ...
                     'comparable.'])
        end
        tempwavfile = isp_tempfile('.wav');
        tempmp3file = isp_tempfile('.mp3');
        wavwrite(wav, options.samplerate, 16, tempwavfile);
        isp_callexecutable('isp_lame', ['--disptime 10 -h --cbr -b ' ...
                            num2str(options.bitrate) ...
                            ' "' tempwavfile '" "' tempmp3file '"']);
        delete(tempwavfile);
        [wav,fs,nbits] = isp_audioread(tempmp3file, 'samplerate', options.samplerate);
        if fs ~= options.samplerate
            fprintf(1, 'Resampling mp3 to original samplerate.\n');
            wav=resample(wav, options.samplerate, fs);
        end
        delete(tempmp3file);
        if nCh == 1
            wav=mean(wav, 2);
        end
    end

    % Mono
    if options.mono
        wav = mean(wav, 2);
    end

    % Change bandwidth
    if ~isinf(options.bandwidth)
        wav = resample(wav, 2*options.bandwidth, options.samplerate);
        wav = resample(wav, options.samplerate, 2*options.bandwidth);
    end

    % Crop to maximum length
    maxsamples = options.maxlength*options.samplerate;
    if maxsamples < size(wav, 1)
        % Restrict length
        startIdx = 1 + floor(.5*(size(wav, 1) - maxsamples));
        wav = wav (startIdx:(startIdx+maxsamples-1), :);
    end
    
    % Add silence
    if options.addsilence ~= 0
        wavLength = size(wav, 1);
        nZeros = wavLength * options.addsilence;
        nChannels = size(wav, 2);
        nZerosBefore = floor(nZeros/2);
        nZerosAfter = nZeros - nZerosBefore;
        noiseAmplitude = 1/32768;
        wav = [noiseAmplitude*randn(nZerosBefore, nChannels)
               wav
               noiseAmplitude*randn(nZerosAfter, nChannels)];
    end

end