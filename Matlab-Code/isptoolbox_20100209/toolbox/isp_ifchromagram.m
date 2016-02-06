function [Y,p,m,S] = isp_ifchromagram(d,sr,fftlen,nbin,f_ctr,f_sd)
%ISP_IFCHROMAGRAM  Compute chromagram from instantaneous frequency
%
% SYNTAX
%   [Y,p,m,S] = isp_ifchromagram(d,sr,fftlen,nbin,f_ctr,f_sd)
%
% DESCRIPTION
%   Calculate the "chromagram" of the sound in d. Use instantaneous
%   frequency to keep only real harmonics. The first time this function
%   is called, it attempts to compile some C code to speed things up. If
%   this fails, a MATLAB version is used instead.
%
% INPUT
%   d:
%     Wave signal.
%   sr:
%     Sampling rate.
%   fftlen:
%     Window length.
%   ffthop:
%     Hop length.
%   nbin:
%     Number of steps to divide the octave into.
%   f_ctr, f_sd:
%     Weight with center frequency f_ctr (in Hz) and gaussian SD f_sd (in
%     octaves).
%
% OUTPUT
%   Y:
%     Chromagram 
%   p:
%     Frequencies of instantaneous frequency gram
%   m:
%     Magnitudes of instantaneous frequency gram
%   S:
%     Complex STFT
%
% SEE ALSO
%   isp_ifgram, isp_ifptrack.
% 
% HISTORY
%   2006-09-26: dpwe@ee.columbia.edu
%   2007-10-11: Heavily optimized by Jesper Højvang Jensen (jhj@es.aau.dk)

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.



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

if nargin < 3;   fftlen = 2048; end
if nargin < 4;   nbin = 12; end
if nargin < 5;   f_ctr = 400; end  % was 1000, but actually I used 400
if nargin < 6;   f_sd = 1; end

A0 = 27.5; % Hz
A440 = 440; % Hz
f_ctr_log = log(f_ctr/A0) / log(2);

fminl = octs2hz(hz2octs(f_ctr)-2*f_sd);
fminu = octs2hz(hz2octs(f_ctr)-f_sd);
fmaxl = octs2hz(hz2octs(f_ctr)+f_sd);
fmaxu = octs2hz(hz2octs(f_ctr)+2*f_sd);

ffthop = fftlen/4;
nchr = nbin;

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


% Calculate spectrogram and IF gram pitch tracks...
[p,m,S]=isp_ifptrack(d,fftlen,sr,fminl,fminu,fmaxl,fmaxu); 

if exist('isp_ifchromagram_helper')==3
    Y = isp_ifchromagram_helper(p, m, nchr);
else
    fprintf(1, 'Couldn''t find mex-file. Using matlab version for %s.\n', mfilename);

    [nbins,ncols] = size(p);

    % chroma-quantized IF sinusoids
    nzp = find(p(:)>0);
    Pocts=p;
    Pocts(nzp) = hz2octs(p(nzp));

    % Figure best tuning alignment
    [hn,hx] = hist(nchr*Pocts(nzp)-round(nchr*Pocts(nzp)),100);
    centsoff = hx(find(hn == max(hn)));
    % The C version of this code is equivalent to the following
    % [hn,hx]= hist(nchr*Pocts(nzp)-round(nchr*Pocts(nzp)),(-49.5:1:50)/100);
    % So, the C version uses 100 bins equally spaced between -0.5 and
    % 0.5, while the hist function uses 100 bins spaced between the
    % minimum and maximum values. The minimum will almost always be
    % really close to -0.5, and the maximum will be really close to 0.5, so the
    % difference is usually negligible. However, since frequencies are
    % rounded to the nearest semitone, even very small rounding errors might
    % cause seemingly large differences.

    % Adjust tunings to align better with chroma
    Pocts(nzp) = Pocts(nzp) - centsoff(1)/nchr;

    % Quantize to chroma bins
    PoctsQ = Pocts;
    PoctsQ(nzp) = round(nchr*Pocts(nzp))/nchr;

    % map IF pitches to chroma bins
    Pmapc = round(nchr*(PoctsQ - floor(PoctsQ)));
    Pmapc(p(:) == 0) = -1; 
    Pmapc(Pmapc(:) == nchr) = 0;

    Y = zeros(nchr,ncols);
    tmp=(0:(nchr-1))';
    for t = 1:ncols;
        Y(:,t)=(tmp(:, ones(1,size(Pmapc,1)))==Pmapc(:,t*ones(nchr,1))')*m(:,t);
    end

end
function octs = hz2octs(freq, A440)
% octs = hz2octs(freq, A440)
% Convert a frequency in Hz into a real number counting 
% the octaves above A0. So hz2octs(440) = 4.0
% Optional A440 specifies the Hz to be treated as middle A (default 440).
% 2006-06-29 dpwe@ee.columbia.edu for fft2chromamx

if nargin < 2;   A440 = 440; end

% A4 = A440 = 440 Hz, so A0 = 440/16 Hz
octs = log(freq./(A440/16))./log(2);

function hz = octs2hz(octs,A440)
% hz = octs2hz(octs,A440)
% Convert a real-number octave 
% into a frequency in Hzfrequency in Hz into a real number counting 
% the octaves above A0. So hz2octs(440) = 4.0.
% Optional A440 specifies the Hz to be treated as middle A (default 440).
% 2006-06-29 dpwe@ee.columbia.edu for fft2chromamx

if nargin < 2;   A440 = 440; end

% A4 = A440 = 440 Hz, so A0 = 440/16 Hz

hz = (A440/16).*(2.^octs);


