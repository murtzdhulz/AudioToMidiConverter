%ISP_TICHROMA  Define time scale invariant chroma-based distance measure.
%
% SYNTAX
%   distancemeasure = isp_tichroma
%
% DESCRIPTION
%   Return a struct that defines a time scale invariant chroma-based
%   distance measure. Functions such as isp_evaluate,
%   isp_extractfeature, isp_computedistance accept this struct as input.
%
% OUTPUT
%   distancemeasure:
%     Struct defining the distance measure.
%
% SEE ALSO
%   isp_evaluate, isp_extractfeature, isp_computedistance, isp_mfccgmmkl.
%
% HISTORY
%   2007: Created by Jesper H. Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function varargout = isp_tichroma(functionName, varargin)
    
    if ~exist('functionName')
        functionName = 'distanceStruct';
    end
    
    varargout{1:nargout} = feval(functionName, varargin{:});
end

function distancemeasure = distanceStruct

    fs = 11025;

    % Define distance measure
    distancemeasure.name = 'Tempoinvariant chroma';
    distancemeasure.computefeature = ...
        ['feature = ' mfilename '(''computefeature'', wav, options);'];
    distancemeasure.computedistance = ...
        ['featureDistance = ' mfilename ...
         '(''computedistance'', feature1, feature2, options);'];
    distancemeasure.computedistancematrix = ...
        ['distancematrix = ' mfilename ...
         '(''computedistancematrix'', features1, features2, options);'];

    % 'usedFunctions' is only needed when using the distributed jobs
    % framework and compiling everything with with mcc.
    distancemeasure.usedFunctions = {mfilename};

    % 'fs' and 'mono' is used by 'isp_extractfeature' to determine the
    % format of the 'wav' variable.
    distancemeasure.samplerate = fs;
    distancemeasure.mono = true;

    % Options passed on to 'computefeature' and
    % 'computedistance[matrix]':
    distancemeasure.options.fs = fs;
    distancemeasure.options.storeChroma = false;
    distancemeasure.options.fs = fs;
    distancemeasure.options.tMax = 18;
    distancemeasure.options.tMin = 1.7;
    distancemeasure.options.nBands = 17;
    distancemeasure.options.maxOffset = 1;
    distancemeasure.options.storeChroma = false;
    distancemeasure.options.acexp = 0.7;
    distancemeasure.options.silent = false;
    distancemeasure.options.variation = 'autocorr';


end

function feature = computefeature(wav, options)
    feature = chromtimeinvftrs(wav, options, 1);
end


function dst = computedistance(feature1, feature2, options)
    dst = computedistancematrix({feature1}, {feature2}, options);
end


function distancematrix = computedistancematrix(features1, features2, options)

    if options.storeChroma
        for n=1:length(features1)
            features1{n} = chromtimeinvftrs(features1{n}, options, 2);
        end
        for n=1:length(features2)
            features2{n} = chromtimeinvftrs(features2{n}, options, 2);
        end
    end

    switch options.variation
      case 'autocorr'
        nTones = size(features1{1}, 1);
        maxOffset = options.maxOffset;

        % Perform the cross-correlation of the two chroma beat ftr matrices
        distancematrix = zeros(numel(features1),numel(features2));
        ffeatures2=zeros(1, numel(features2));
        for n=1:numel(features2)
            tmp = [zeros(nTones, maxOffset) ...
                             features2{n} ...
                             zeros(nTones, maxOffset)];
            ffeatures2(1:numel(tmp),n) = tmp(:);
        end
        tmp=zeros(2*maxOffset+1, size(ffeatures2,1));
        for iFeature=1:numel(features1)
            for n=0:2*maxOffset
                tmp(n+1,:) = [zeros(nTones*n,1)
                               features1{iFeature}(:)
                               zeros(nTones*(2*maxOffset-n),1)]';
            end
            distancematrix(iFeature, :) = 2-2*max(tmp*ffeatures2, [], 1);
        end
        distancematrix(distancematrix(:) < 0) = 0;

      otherwise
        nTones = size(features1{1}, 1);
        maxOffset = options.maxOffset;
        L = max(size(features1{1}, 2), size(features2{1}, 2)) + maxOffset;

        % Avoid unfortunate fft sizes
        if 27 <= L && L < 32
            L = 32;
        end
        if L == 37, L = 38; end
        if L == 41, L = 42; end
        if L == 43, L = 44; end
        if L == 47, L = 48; end


        ffeatures1 = cell(numel(features1), 1);
        ffeatures2 = cell(numel(features2), 1);
        for n=1:numel(features1)
            ffeatures1{n} = fft2(features1{n}, nTones, L);
        end
        for n=1:numel(features2)
            ffeatures2{n} = fft2(features2{n}, nTones, L);
        end

        distancematrix = zeros(numel(features1),numel(features2));
        validOffsets = [1:maxOffset+1 size(ffeatures1{1}, 2)-(0:maxOffset-1)];

        % Perform the cross-correlation of the two chroma beat ftr matrices
        for iFeature=1:numel(features1)
            if ~options.silent
                fprintf(1, 'Computing distance matrix row %d of %d.\n', ...
                        iFeature, numel(features1));
            end
            for jFeature=1:numel(features2)
                distances = 1 - ifftn(ffeatures1{iFeature}.* ...
                                      conj(ffeatures2{jFeature}), 'symmetric');
                distancematrix(iFeature,jFeature) = min(min(distances(:, validOffsets)));
            end
        end

        % Rounding errors might cause negative distances. Truncate to zero.
        distancematrix(distancematrix(:) < 0) = 0;
    end    
end









function [F,bts] = chromtimeinvftrs(d,options, step, f_ctr,f_sd)
%
% This function is derived from Dan Ellis' chrombeatftrs.m and is thus
% released under the GPL.
%
% Help text from the original function:
%
% [F,bts] = chrombeatftrs(D,SR,F_CTR,F_SD,TYPE)
%    F returns a feature vector of beat-level chroma features (12
%    rows x n time step columns). bts returns the times of all the
%    beats.  
%    New version separates out chroma calculation
%    TYPE selects chroma calculation type; 1 (default) uses IF; 
%    2 uses all FFT bins, 3 uses only local peaks (a bit like Emilia).
% 2006-07-14 dpwe@ee.columbia.edu

%   Copyright (c) 2006 Columbia University and 2007 Jesper Højvang Jensen.
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

    if nargin < 4; f_ctr = 1000; end
    if nargin < 5; f_sd = 1; end


    sr = options.fs;
    tMax = options.tMax;
    tMin = options.tMin;
    nBands = options.nBands;

    % Calculate frame-rate chromagram
    fftlen = 2 ^ (round(log(sr*(2048/22050))/log(2)));
    nbin = 12;

    if ~options.storeChroma || (options.storeChroma && step==1)

        Y = isp_ifchromagram(d,sr,fftlen,nbin,f_ctr,f_sd);
        
        if size(Y, 2) == 0
            warning('Audio file too short.')
            Y = zeros(12, 1);
        end

        if options.storeChroma
            F=Y;
            return
        end
    else
        Y=d;
    end
    

    ffthop = fftlen/4;
    fs = sr / ffthop; % Chromagram samplerate
    switch options.variation
        case {'logfreq', 'sqrtfreq'}
          switch options.variation
            case 'logfreq'
              delta = 0.01;
              freqChroma = abs(fft(log((Y+delta)/delta), [], 2)).^2;
            case 'sqrtfreq'
              freqChroma = abs(fft(sqrt(Y), [], 2)).^2;
            otherwise
              error('This should never happen :-(')
          end

          flen = size(freqChroma, 2);
          
          boundaries = fliplr(1./logspace(log10(tMin), log10(tMax), nBands+2));
          
          lowerFrequency = boundaries(1:end-2);
          centerFrequency = boundaries(2:end-1);
          upperFrequency = boundaries(3:end);
          
          frequencies = fs / flen * (0:flen-1);
          
          BandMatrix=zeros(nBands, length(frequencies));
          for iBand=1:nBands
              betweenLowerAndCenter = frequencies > lowerFrequency(iBand) ...
                  & frequencies <= centerFrequency(iBand);
              betweenCenterAndUpper = frequencies < upperFrequency(iBand) ...
                  & frequencies > centerFrequency(iBand);
              BandMatrix(iBand, :) = betweenLowerAndCenter ...
                  .* (frequencies - lowerFrequency(iBand)) ...
                  / (centerFrequency(iBand)-lowerFrequency(iBand));
              BandMatrix(iBand, :) = BandMatrix(iBand, :) + ...
                  betweenCenterAndUpper .* (frequencies - upperFrequency(iBand)) ...
                  / (centerFrequency(iBand)-upperFrequency(iBand));
          end
          
          F = freqChroma * BandMatrix';
      
          
      case 'autocorr'
        Y=Y.^options.acexp;
        
        %Y = max(0,diff(Y, 1, 2));

        Y = filter([1 -1], [1 -.99], Y, [], 2);
        minFftLength = tMax*fs + size(Y,2);
        fftLength = 2^ceil(log2(minFftLength));
        
        tmp = fftn(Y, [size(Y,1) fftLength]);
        autocorr = ifftn(tmp.*conj(tmp));
        boundaries = logspace(log10(tMin), log10(tMax), nBands+2);

        startTime = boundaries(1:end-2);
        centerTime = boundaries(2:end-1);
        endTime = boundaries(3:end);
        lags = (0:size(autocorr,2)-1)/fs;
        autocorr(:, lags>endTime(end)) = [];
        lags(lags>endTime(end)) = [];
        BandMatrix=zeros(nBands, length(lags));
        for iBand=1:nBands
            betweenLowerAndCenter = lags > startTime(iBand) ...
                & lags <= centerTime(iBand);
            betweenCenterAndUpper = lags < endTime(iBand) ...
                & lags > centerTime(iBand);
            BandMatrix(iBand, :) = betweenLowerAndCenter ...
                .* (lags - startTime(iBand)) ...
                / (centerTime(iBand)-startTime(iBand));
            BandMatrix(iBand, :) = BandMatrix(iBand, :) + ...
                betweenCenterAndUpper .* (lags - endTime(iBand)) ...
                  / (centerTime(iBand)-endTime(iBand));
        end
        BandMatrix = diag(1./sum(BandMatrix, 2)') * BandMatrix;

        F = autocorr * BandMatrix';

      otherwise
        error('Invalid value for options.variation.')
    end

    p = norm(F(:));
    F = F/(p + (p==0));
end




