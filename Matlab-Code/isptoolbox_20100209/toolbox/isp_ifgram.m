%ISP_IFGRAM  Instantaneous frequency by phase derivative.
%
% SYNTAX
%   [F,D] = ifgram(X, N, W, H, SR, maxbin)
%
% DESCRIPTION
%    Compute the instantaneous frequency (as a proportion of the sampling
%    rate) obtained as the time-derivative of the phase of the complex
%    spectrum as described by Toshihiro Abe et al in ICASSP'95,
%    Eurospeech'97. Calculates regular STFT as side effect.
%
% INPUT
%   X:
%     Wave signal.
%   N:
%     FFT length.
%   W:
%     Window length.
%   H:
%     Step length.
%   SR:
%     Sampling rate.
%   maxbin:
%     The index of the maximum bin needed. If specified, unnecessary
%     computations are skipped.
%
% OUTPUT
%   F:
%     Instantaneous frequency spectrogram.
%   D:
%     Short time Fourier transform spectrogram.
%
% SEE ALSO
%   isp_ifchromagram, isp_ifgram, isp_ifptrack.
%
% HISTORY
%   after 1998may02:  dpwe@icsi.berkeley.edu
%   2001-03-05:  dpwe@ee.columbia.edu  revised version
%   2001-12-13:  dpwe@ee.columbia.edu  Fixed to work when N != W
%   2007-2008:  Optimized by Jesper Højvang Jensen
%   2008-08-14:  Memory leak fixed by Dan Ellis.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.



function [F,D] = isp_ifgram(X, N, W, H, SR, maxbin)
% Dan's original copyright:
%
%   Copyright (c) 2006 Columbia University.
% 
%   This file is part of LabROSA-coversongID
% 
%   LabROSA-coversongID is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License version 2 as
%   published by the Free Software Foundation.
% 
%   LabROSA-coversongID is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with LabROSA-coversongID; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
%   02110-1301 USA
% 
%   See the file "COPYING" for the text of the license.

if nargin < 2;  N = 256; end
if nargin < 3;  W = N;   end
if nargin < 4;  H = W/2; end
if nargin < 5;  SR = 1;  end
if nargin < 6;  maxbin = 1+N/2; end


% Compile mex file if necessary
mexname = [mfilename '_helper'];
if ~(exist(mexname)==3)
    if exist('OCTAVE_VERSION')
        npname = fullfile(isp_toolboxpath, mexname);
        mkoctfile([npname '.c'], '--mex', '--output', [npname '.' mexext])
    else
        cfile=which([mexname '.c']);
        try
            mex(cfile, '-outdir', fileparts(cfile))
        catch
            warning('Unable to compile mex file')
        end
    end
end

if ~isvector(X)
    error('X must be a vector')
end

if exist('isp_ifgram_helper')==3
    [F,D] = isp_ifgram_helper(X, N, W, H, SR, maxbin);
else
    fprintf(1, 'Couldn''t find mex-file. Using matlab version for %s.\n', mfilename);
    Flen = maxbin;
    s = length(X);
    % Make sure it's a single row
    if size(X,1) > 1
        X = X';
    end
    %win = [0,hanning(W-1)'];
    win = 0.5*(1-cos([0:(W-1)]/W*2*pi));

    % Window for discrete differentiation
    T = W/SR;
    dwin = -pi / T * sin([0:(W-1)]/W*2*pi);

    % sum(win) takes out integration due to window, 2 compensates for neg frq
    norm = 2/sum(win);

    % How many complete windows?
    nhops = 1 + floor((s - W)/H);

    F = zeros(Flen, nhops);
    D = zeros(Flen, nhops);

    nmw1 = floor( (N-W)/2 );
    %nmw2 = N-W - nmw1;

    %ww = 2*pi*[0:(N-1)]*SR/N;
    ww = 2*pi*[0:(Flen-1)]*SR/N;

    wu = zeros(1, N);
    du = zeros(1, N);

    for h = 1:nhops
        u = X((1+(h-1)*H):(W+(h-1)*H));

        % Pad or truncate samples if N != W
        % Apply windows now, while the length is right
        if N >= W
            wu(nmw1+1:nmw1+W) = win.*u;
            du(nmw1+1:nmw1+W) = dwin.*u;
        elseif N < W
            wu = win(1-nmw1:N-nmw1).*u(1-nmw1:N-nmw1);
            du = dwin(1-nmw1:N-nmw1).*u(1-nmw1:N-nmw1);
        end

        % FFTs of straight samples plus differential-weighted ones
        % Replaced call to fftshift with inline version. Jesper Højvang Jensen, Aug 2007
        % t1 = fft(fftshift(du));
        % t2 = fft(fftshift(wu));
        split = ceil(length(du)/2) + 1;
        t1 = fft(du([split:end 1:split-1]));
        t2 = fft(wu([split:end 1:split-1]));
        t1 = t1(1:Flen);
        t2 = t2(1:Flen);


        % Scale down to factor out length & window effects
        D(:,h) = t2'*norm;

        % Calculate instantaneous frequency from phase of differential spectrum
        t = t1 + j*(ww.*t2);
        a = real(t2);
        b = imag(t2);
        da = real(t);
        db = imag(t);
        instf = (1/(2*pi))*(a.*db - b.*da)./((a.*a + b.*b)+(t2==0));
        % 1/2pi converts rad/s into cycles/s
        % sampling rate already factored in as constant in dwin & ww
        % so result is in Hz
        
        F(:,h) = instf';

    end;

end