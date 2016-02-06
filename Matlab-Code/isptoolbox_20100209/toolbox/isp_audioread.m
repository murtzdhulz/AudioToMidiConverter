function [Y,FS,NBITS,OPTS]=isp_audioread(filename, varargin)
%ISP_AUDIOREAD  Read audio file.
%
% SYNTAX
%   [y,fs,nbits,opts]=isp_audioread(filename, options, 'field1', value1, ...)
%
% DESCRIPTION
%   ISP_READAUDIO calls mplayer to convert media files supported by mplayer
%   to .wav, which is then returned as the functions output. The mplayer
%   binary must be obtained separately and be placed either in the ISP
%   toolbox directory or in the path.
%   
% INPUT 
%   filename:
%     Name of audio file.
%   options, field/value pairs:
%     The following parameters can be set as field names in the
%     'options' struct, or be specified as field/value pairs:
%
%     nChannels:
%       Number of output channels. If it is empty, use the number of
%       channels in the input (Default: empty).
%     samplerate:
%       Sampling frequency of Y. If it is empty, use the sample
%       rate of the input file. (Default: empty).
%     maxSamples:
%       The maximum number of samples to read. Only intended to avoid
%       "Out of memory" erros when handling very large audio
%       files. (Default: inf)
%  
% OUTPUT
%   y:
%     Wave signal.
%   fs:
%     Sampling rate of Y.
%   nbits:
%     Number of bits.
%   opts:
%     Struct with additional information as returned by WAVREAD (not
%     available with Octave at the time of writing).
%
% EXAMPLE
%   Read the file 'song.mp3':
%     [y,fs] = isp_audioread('song.mp3');
%   The following all read the file 'song.mp3' as 16 kHz mono:
%     y1 = isp_audioread('song.mp3', 'samplerate', 16000, 'nChannels', 1);
%     
%     opts.samplerate=16000;
%     y2 = isp_audioread('song.mp3', opts, 'nChannels', 1);
%     
%     opts.nChannels=1; opts.samplerate=16000;
%     y3 = isp_audioread('song.mp3', opts);
%
% SEE ALSO
%   wavread
%
% HISTORY
%   Created by Jesper H. Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.



% Switches for MPlayer configure:
% Misc.:
% --enable-runtime-cpudetection
% Audio:
% --disable-alsa --disable-ossaudio --disable-arts --disable-esd --disable-polyp --disable-jack --disable-openal --disable-nas --disable-sgiaudio --disable-sunaudio --disable-win32waveout --disable-select
% Video:
% --disable-vidix-internal --disable-vidix-external --disable-gl --disable-dga --disable-vesa --disable-svga --disable-sdl --disable-aa --disable-caca --disable-ggi --disable-ggiwmh --disable-directx --disable-dxr2 --disable-dxr3 --disable-ivtv --disable-dvb --disable-mga --disable-xmga --disable-xv --disable-xvmc --disable-vm --disable-xinerama --disable-x11 --disable-xshape --disable-fbdev --disable-mlib --disable-3dfx --disable-tdfxfb --disable-s3fb --disable-directfb --disable-zr --disable-bl --disable-tdfxvid --disable-tga --disable-pnm --disable-md5sum
% Optional features. Had my doubts about disabling smb ...
% --disable-mencoder --disable-termcap --disable-termios --disable-iconv --disable-langinfo --disable-lirc --disable-lircc --disable-vm --disable-xf86keysym --disable-radio-v4l2 --disable-tv --disable-tv-v4l1 --disable-tv-v4l2 --disable-tv-bsdbt848 --disable-pvr --disable-rtc --disable-winsock2 --disable-live --disable-dvdnav --disable-dvdread --disable-mpdvdkit --disable-cdparanoia --disable-bitmap-font --disable-freetype --disable-fontconfig --disable-unrarlib --disable-sortsub --disable-fribidi --disable-enca --disable-macosx --disable-maemo --disable-macosx-bundle --disable-vstream --disable-pthreads --disable-ass --disable-smb
%
% All switches were thrown after configure before running make:
% ./configure  --enable-runtime-cpudetection  --disable-alsa --disable-ossaudio --disable-arts --disable-esd --disable-polyp --disable-jack --disable-openal --disable-nas --disable-sgiaudio --disable-sunaudio --disable-win32waveout --disable-select  --disable-vidix-internal --disable-vidix-external --disable-gl --disable-dga --disable-vesa --disable-svga --disable-sdl --disable-aa --disable-caca --disable-ggi --disable-ggiwmh --disable-directx --disable-dxr2 --disable-dxr3 --disable-ivtv --disable-dvb --disable-mga --disable-xmga --disable-xv --disable-xvmc --disable-vm --disable-xinerama --disable-x11 --disable-xshape --disable-fbdev --disable-mlib --disable-3dfx --disable-tdfxfb --disable-s3fb --disable-directfb --disable-zr --disable-bl --disable-tdfxvid --disable-tga --disable-pnm --disable-md5sum  --disable-mencoder --disable-termcap --disable-termios --disable-iconv --disable-langinfo --disable-lirc --disable-lircc --disable-vm --disable-xf86keysym --disable-radio-v4l2 --disable-tv --disable-tv-v4l1 --disable-tv-v4l2 --disable-tv-bsdbt848 --disable-pvr --disable-rtc --disable-winsock2 --disable-live --disable-dvdnav --disable-dvdread --disable-mpdvdkit --disable-cdparanoia --disable-bitmap-font --disable-freetype --disable-fontconfig --disable-unrarlib --disable-sortsub --disable-fribidi --disable-enca --disable-macosx --disable-maemo --disable-macosx-bundle --disable-vstream --disable-pthreads --disable-ass --disable-smb && make clean && make
% Furthermore, libmad.so.0 had to be copied to the isptoolbox directory
% as well.    
    

    options.nChannels = [];
    options.samplerate = [];
    options.maxSamples = inf;
    options = isp_interpretarguments(options, varargin{:});
  
    % Give an error if a regular file is specified (i.e., no prefix such
    % as http:) and it doesn't exist
    if isempty(regexp(filename, '^[a-z]*:')) && ~exist(filename, 'file')
        error('File does not exist');
    end

    if any(exist('OCTAVE_VERSION') == [5 102])
        isOctave = true;
    else
        isOctave = false;
    end
  
    % Decode audio file
    tempwavfile = isp_tempfile('.wav');

    executedir='';
    if strcmp(computer, 'PCWIN') && length(tempwavfile) > 1 && tempwavfile(2)==':'
        % The colon from Windows drive letters is confused with yet a
        % pcm option
        executedir=tempwavfile(1:2);
        if length(filename)<2 || filename(2)~=':'
            c=pwd;
            filename=[c(1:2) filename];
        end
        mplayerOpt = ['-noconsolecontrols -vc null -vo null -quiet -ao pcm:file=' tempwavfile(3:end) ':waveheader:fast'];
        pre='';
    else
        mplayerOpt = ['-noconsolecontrols -vc null -vo null -quiet -ao pcm:file=' tempwavfile ':waveheader:fast'];
        % Need this one for libmad.so.0
        pre='LD_LIBRARY_PATH="{}" ';
    end        

    if ~isempty(options.samplerate)
        mplayerOpt = [mplayerOpt ' -af resample=' num2str(options.samplerate) ':1:1'];
    end

    isp_callexecutable('isp_mplayer', [ mplayerOpt ' "' filename '"'], ...
                       executedir, pre);

    if isOctave
        [Y,FS,NBITS] = wavread(tempwavfile);
    else
        % An extra check to ensure the end of the file was not cut off due
        % to lack of disk space. Matlabs wavread only complains about this
        % if the 'size' mode is used, not otherwise.
        Ysize = wavread(tempwavfile, 'size');
        if Ysize(1) > options.maxSamples
            warning('Truncating audio file')
            [Y,FS,NBITS,OPTS] = wavread(tempwavfile, options.maxSamples);
            if size(Y,1) ~= options.maxSamples
                error('Error reading decoded MP3 file.')
            end
        else
            [Y,FS,NBITS,OPTS] = wavread(tempwavfile);
            if ~isequal(Ysize, size(Y))
                error('Error reading decoded MP3 file.')
            end
        end
    end
    delete(tempwavfile);

    if ~isempty(options.nChannels)
        if options.nChannels == 1
            Y = mean (Y, 2);
        else
            if options.nChannels ~= size(Y,2)
                warning(['Sorry, unable to provide the requested number ' ...
                         'of channels'])
            end
        end
    end
end