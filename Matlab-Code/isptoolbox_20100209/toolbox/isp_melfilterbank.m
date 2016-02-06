%ISP_MELFILTERBANK  Compute the mel scaled filterbank.
%
% SYNTAX
%   [H,options,fk] = isp_melfilterbank(options ...)
%
% DESCRIPTION
%   Compute the mel scaled filterbank.
%
% INPUT
%   fs:
%     Sampling frequency (default: 44100).
%   tw:
%     Window size in seconds (default: 0.02).
%   fl:
%     Lowest frequency of the filter bank (default: 0).
%   fh:
%     Highest frequency of the filter bank (default: 11025).
%   Nfilters:
%     Number of filters in the filterbank (default: 30).
%
% OUTPUT
%   H:
%     Mel-scaled filterbank.
%   options:
%     Input options supplemented with default values.
%   fk:
%     Spectrogram frequencies.
%
% SEE ALSO
%   isp_mfccsig.
%
% HISTORY
%   Created by Sigurdur Sigurdsson, ISP, IMM, DTU.
%   Modified for toolbox by Jesper H. Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


% Set default values if necessary

function [H,options,fk] = isp_melfilterbank(varargin)

options = isp_interpretarguments(struct(...
    'Nfilters', 30, ...
    'fs', 44100, ...
    'fh', 11025, ...
    'fl', 0, ...
    'tw', 0.02), varargin{:})


Nfilters = options.Nfilters;
fs = options.fs;
fh = options.fh;
fl = options.fl;
tw = options.tw;

% Compute the number of fft points
Nfft = round(tw*fs);

% Compute the center frequencies for the Mel filters
phi_h = 2595*log10(fh/700+1);
phi_l = 2595*log10(fl/700+1);
d_phi = (phi_h-phi_l)/(Nfilters+1);
phi_c = d_phi*(0:(Nfilters+1));
fc = 700*(10.^((phi_c/2595))-1);

% Compute the frequencies of the spectrogram
fk = ((0:(Nfft-1))/Nfft)*fs;
fk = fk(1:(floor(Nfft/2)+1));

% Compute the Mel filter bank with scaled energy
H = zeros(Nfilters,(floor(Nfft/2)+1));
for m = 2:(Nfilters+1)
  idx1 = fk>fc(m-1) & fk< fc(m);
  idx2 = fk>=fc(m) & fk<= fc(m+1);
  H(m-1,idx1) = (fk(idx1)-fc(m-1))/(fc(m)-fc(m-1));
  H(m-1,idx2) = (fk(idx2)-fc(m+1))/(fc(m)-fc(m+1));
  H(m-1,:) = H(m-1,:)/(2*(fc(m+1)-fc(m-1)));
end

% If no output, then plot filter bank
if nargout == 0
  plot(fk,H,'.-')
  xlabel('Frequency (Hz)')
  ylabel('Amplitude')
end