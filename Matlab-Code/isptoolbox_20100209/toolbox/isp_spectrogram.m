function [X,options,f,t] = isp_spectrogram(x,options)
%ISP_SPECTROGRAM  Compute the spectrogram of a signal.
%
% SYNTAX:
%   [X,options,f,t] = isp_spectrogram(x,options)
%
% DESCRIPTION
%   The function isp_spectrogram computes the spectrogram of a signal.
%   The spectrogram is computed using the asinh of the absolute discrete
%   Fourier transform, where a Hamming window is applied.
%
% INPUT
%   x:
%     Signal vector.
%   options ...:
%     Structs or field/value pairs with any of the following fields:
%     fs:
%       Sampling frequency (default 44100).
%     tw:
%       Window size in seconds (default 0.02).
%     to:
%       Overlap in seconds (default 0.01)        .
%
% OUTPUT
%   X:
%     Spectrogram for signal x
%   options:
%     Options used to compute the spectrogram
%   f:
%     Frequency axis for spectrogram
%   t:
%     Time axis for spectrogram
%
% SEE ALSO
%   isp_mfccsig.
%
% HISTORY
%   Created by Sigurdur Sigurdsson, ISP, IMM, DTU.
%   Modified to suit toolbox by Jesper H. Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


% Set default values if necessary
if nargin == 1, options.mfcc_spectrogram = struct; end
if ~isfield(options,'spectrogram'), options.spectrogram = struct; end
if ~isfield(options.spectrogram,'tw'), options.tw=0.02; end
if ~isfield(options.spectrogram,'fs'), options.fs=44100; end
if ~isfield(options.spectrogram,'to'), options.to=0.01; end

% Get option values
tw = options.tw;
to = options.to;
fs = options.fs;

% Check for one dimensional data
[N,D] = size(x);
if min(N,D)~=1
  error('Only one dimensional signal may be applied.')
end

% Make a column vector of the data
x = x(:);

% Compute the length of window in samples
N = length(x);
Nfft = round(tw*fs);

% Compute the start and stop indexes for all clips
Nover = Nfft-round(to*fs);
idx_start = 1:Nover:N;
idx_stop = (Nfft):Nover:N;
K = min(length(idx_start),length(idx_stop));
idx_start = idx_start(1:K);
idx_stop = idx_stop(1:K);

% Organize the signal into clips and Hamming window all clips
w = hamming(Nfft);
x1 = zeros(Nfft,K);
for k = 1:K
  x1(:,k) = x(idx_start(k):idx_stop(k)).*w;
end

% Compute the DFT of the clips
X = fft(x1);

% Compute the absolute value and scale with asinh
%X = asinh(abs(X(1:(floor(Nfft/2)+1),:)));
X = abs(X(1:(floor(Nfft/2)+1),:));
X = asinh(X);

% Compute frequency axis
f = ((0:(Nfft-1))/Nfft)*fs;
f = f(1:(floor(Nfft/2)+1));

% Compute time axis
t = (idx_stop-idx_start)/(2*fs);

% If no output, then plot the spectrogram
if nargout == 0
  imagesc(t,f,X)
  set(gca,'Ydir','Normal')
  colorbar
  xlabel('Time (sec)')
  ylabel('Frequency (Hz)')
end