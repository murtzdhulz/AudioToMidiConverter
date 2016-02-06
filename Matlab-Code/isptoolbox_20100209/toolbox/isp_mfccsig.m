%ISP_MFCCSIG  Compute the MFCC of a signal.
%
% SYNTAX
%   [mfcc, options] = isp_mfcc(x,options)
%
% DESCRIPTION
%   Compute the MFCC of a signal.
% 
% INPUT
%   x:
%     Wave signal.
%   options:
%     Struct with any of the following fields:
%     fs:
%       Sampling frequency (default 44100).
%     tw:
%       Window size in seconds (default 0.02).
%     to:
%       Window overlap in seconds (default 0.01).
%     fl:
%       Lowest frequency of the filter bank (default 0).
%     fh:
%       Highest frequency of the filter bank (default 11025).
%     Nfilters:
%       Number of filters in the filterbank (default 30).
%     Nmfcc:
%       Number of returned MFCCs (default 10).
%
% OUTPUT
%   mfcc:
%     Spectrogram for signal x
%   options:
%     Input options supplemented with default values.
%
% HISTORY
%   Created by Sigurdur Sigurdsson, ISP, IMM, DTU.
%   Modified to suit toolbox by Jesper H. Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [mfcc,options] = isp_mfccsig(x,options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1, options = struct; end

if ~isfield(options,'Nfilters'), 
  options.Nfilters = 30; 
  options.melfilterbank.Nfilters = 30; 
else
  options.melfilterbank.Nfilters = options.Nfilters;   
end

if ~isfield(options,'Nmfcc'), 
  options.Nmfcc = 10; 
  options.dct.Nmfcc = 10; 
else
  options.dct.Nmfcc = options.Nmfcc;   
end

if ~isfield(options,'fs'),
  options.fs = 44100;
  options.spectrogram.fs = 44100;
  options.melfilterbank.fs = 44100; 
else
  options.spectrogram.fs = options.fs;
  options.melfilterbank.fs = options.fs;   
end

if ~isfield(options,'fh'), 
  options.fh = 11025; 
  options.melfilterbank.fh = 11025;
else
  options.melfilterbank.fh = options.fh;
end

if ~isfield(options,'fl'),
  options.fl = 0;
  options.melfilterbank.fl = 0; 
else
  options.melfilterbank.fl = options.fl; 
end

if ~isfield(options,'tw'),
  options.tw = 0.02;
  options.spectrogram.tw = 0.02;
  options.melfilterbank.tw = 0.02;
else
  options.spectrogram.tw = options.tw;
  options.melfilterbank.tw = options.tw;  
end

if ~isfield(options,'to'),
  options.to = 0.01; 
  options.spectrogram.to = 0.01;
else
  options.spectrogram.to = options.to;  
end

if ~isfield(options.dct,'Nmfcc'), 
  options.Nmfcc = 30; 
  options.dct.Nmfcc = 30;
else
  options.dct.Nmfcc = options.Nmfcc;  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code starts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.fs/2<options.fh
  error('The highest frequency of the Mel filter bank (fh) is higher than half the sampling frequency (fs/2).')
end

% Compute the spectrogram
[X,options] = isp_spectrogram(x,options.spectrogram);

% Generate the Mel scaled filter bank
[H,options] = isp_melfilterbank(options.melfilterbank);

% Compute the Mel scaled spectrogram
Xmel = H*X;

% Obtain the MFCCs by taking applying the DCT of the Mel scaled spectrogram
[mfcc,options] = isp_dct(Xmel,options);





function [mfcc,options] = isp_dct(Xmel,options)
% The function isp_dct computes real the discrete cosine tranform (DCT).
%
% @examples
% mfcc = isp_dct(Xmel)
% [mfcc,options] = isp_dct(Xmel,options)
%
% @options
%  dct.Nmfcc : Number of retained MFCCs (default 10)
%
% @function [mfcc,options] = isp_dct(Xmel,options);
%
% @input Xmel - Mel scaled spectrogram (or any signal)
%        options - options for DCT compution
%
% @output mfcc - Mel frequency cepstral coefficients x
%         options - options for DCT computation
%
% @author Sigurdur Sigurdsson, ISP, IMM, DTU.
% @version 1.0

[Nmfcc,N] = size(Xmel);

if nargin == 0, options.dct = struct; end
if ~isfield(options,'dct'), options.dct = struct; end
if ~isfield(options.dct,'Nmfcc'), options.dct.Nmfcc = 10; end

if Nmfcc < options.dct.Nmfcc
  error('The number of MFCC is larger then the number of Mel filterbanks.')
end
  
Nmfcc = options.dct.Nmfcc;

mfcc = (fft(Xmel,2*(Nmfcc)));
mfcc = (mfcc(1:Nmfcc,:));
scale = (exp(-i*(0:Nmfcc-1)*pi/(2*Nmfcc))).';
mfcc = real(mfcc.*scale(:,ones(1,N)));

% A elementwise computation of the DCT, corresponding to the above code
% n = (0:(Nmfcc-1))';
% for m = 1:size(Xmel,2)
%   for k = 0:(Nmfcc-1)
%     mfcc(k+1,m) = sum(Xmel(:,m).*cos(pi/Nmfcc*(n+0.5)*k));
%   end
% end
