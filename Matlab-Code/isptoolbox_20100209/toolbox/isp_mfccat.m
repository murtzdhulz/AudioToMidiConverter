%ISP_MFCCAT  Reimplementation of MFCC from the Auditory Toolbox
%
% SYNTAX
%   [mfcc, options]=isp_mymfcc3(wav, options)
%
% DESCRIPTION
%   
%   As mymfcc3.m, but uses log(|x(f)| + k) instead of log(|x(f)|) and
%   optionally ignores zero-frames.
%
% INPUT
%   wav:
%     Wave signal.
%   options ...:
%     Structs or field/value pairs with any of the following fields:
%     nLinearBands:
%       The number of linear bands. Default: 13.
%     nLogarithmicBands:
%       The number of logarithmic bands above the linear bands. Default: 27.
%     linearBandwidth:
%       Bandwith of the first 'nLinearBands' linear bands. Default: 66.7.
%     logarithmicScale:
%       Scaling factor of the logarithmic bands, i.e., band N+1 has
%       bandwidth 'logarithmicScale' times the bandwidth of band N.
%       Default: 1.0711703.
%     lowestFrequency:
%       Starting frequency of the first band. Default: 133.
%     nDctCoefficients:
%       Number of DCT coefficients to keep. Default: 13.
%     samplerate:
%       Sample rate of input signal.
%     windowSize:
%       Window size. Default: round(samplerate * 512 / 22050).
%     stepSize:
%       Hop size. Default: round(samplerate * 256 / 22050).
%     fftSize:
%       FFT size. Default: max(512, 2*windowSize).
%     spectrumType:
%       Spectrum type. Possibilities: 'fft', 'mvdr', 'lpc', 'warpedMvdr',
%       'warpedLpc', specifying the use of the fast Fourier transform,
%       minimum variance distortionless response (Capon), linear
%       prediction or warped LPC/MVDR for spectral estimation.
%     preemphasis:
%       Boolean specifying whether to apply a preemphasis
%       filter. Default: True.
%     hammingWindow:
%       Boolean specifying whether to use a rectangular window or a
%       Hamming window. Default: True.
%     logAmplitude:
%       Boolean specifying whether to take the logarithm of the power
%       spectrum. Default: True.
%     minAmplitude:
%       Add a constant k given by k = minAmplitude*sqrt(mean(wav.^2)) when
%       taking the logarithm, i.e., log(x + k) instead of log(x).
%     ignoreSilentFrames:
%       Remove frames with silence.
%
% OUTPUT
%   mfcc:
%     The computed MFCCs.
%   options:
%     Input options supplemented by default values.
%
% HISTORY
%   By Jesper H. Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [mfcc, options, debug]=isp_mfccat(wav, varargin)

    
    function warpedF=warpFrequency(f)
        normalized = 2*pi*f/samplingFrequency;
        warpedNorm = normalized + 2*atan(warpFactor*sin(normalized)./(1-warpFactor*cos(normalized)));
        warpedF = warpedNorm * samplingFrequency/(2*pi);
    end

    function [R, maxamp, maxidx]=myxcorr(data, maxlag)
        f=fft(data, fftSize);
        R=ifft(f.*conj(f));
        if useWindow
            % Unbiased estimate
            N=length(data);
            R=R(1:maxlag+1)./((N:-1:N-maxlag)');
        else
            % Biased estimate
            R=R(1:maxlag+1)/length(data);
        end
        if nargout > 1
            [maxamp, maxidx]=max(f(first:last).*conj(f(first:last)));
            maxamp=sqrt(maxamp/fftSize);
        end
    end

    if ~exist('options', 'var') || ~isfield(options, 'mfcc')
        options = struct;
    end
    
    options = isp_interpretarguments(struct(...
                       'nLinearBands', 13, ...
                       'nLogarithmicBands', 27, ...
                       'linearBandwidth', 66.66666666, ...
                       'logarithmicScale', 1.0711703, ...
                       'lowestFrequency', 133.3333, ...
                       'nDctCoefficients', 13, ...
                       'samplerate', 0, ...
                       'windowSize', nan, ...
                       'stepSize', nan, ...
                       'fftSize', nan, ...
                       'spectrumType', 'fft', ...
                       'preemphasis', true, ...
                       'hammingWindow', true, ...
                       'logAmplitude', true, ...
                       'minAmplitude', 0.001, ...
                       'ignoreSilentFrames', true), ...
                                     varargin{:});

    samplingFrequency = options.samplerate;

    if samplingFrequency == 0
        error('Sampling frequency not specified.');
    end

    if isnan(options.windowSize)
        options.windowSize=round(samplingFrequency * 512 / 22050);
    end

    if isnan(options.stepSize)
        options.stepSize =round(samplingFrequency * 256 / 22050);
    end

    nLinearBands = options.nLinearBands;
    nLogarithmicBands = options.nLogarithmicBands;
    linearBandwidth = options.linearBandwidth;
    logarithmicScale = options.logarithmicScale;
    lowestFrequency = options.lowestFrequency;
    nDctCoefficients = options.nDctCoefficients;
    windowSize = options.windowSize;
    stepsize = options.stepSize;
    spectrumMethod = options.spectrumType;
    minConst = options.minAmplitude*sqrt(mean(wav.^2));
    ignoreSilent = options.ignoreSilentFrames;

    if isnan(options.fftSize)
        % If fftSize is much less than 512, the frequency domain
        % resolution will (probably) be too low (at least at 22 kHz - at
        % 8 kHz, things might need to be adapted)
        options.fftSize = max(512, 2*windowSize);
    end
    fftSize = options.fftSize;

    if ~strcmp(spectrumMethod, 'fft')
        mvdrOrder=options.modelOrder;
    end

    preemphasis = options.preemphasis;
    useWindow = options.hammingWindow;
    takeLog = options.logAmplitude;
    nBands = nLinearBands + nLogarithmicBands;
    nCols=floor((length(wav)-(windowSize-stepsize))/stepsize);

    if isempty(regexp(spectrumMethod, '^warped'))
        warped=false;
    else
        warped=true;
    end

    
    %
    % Specify bands
    %    
    
    lowerFrequency(1:nLinearBands) = lowestFrequency + (0:nLinearBands-1)*linearBandwidth;
    lowerFrequency(nLinearBands+1:nBands) = lowerFrequency(nLinearBands) * logarithmicScale.^(1:nLogarithmicBands);

    centerFrequency = [lowerFrequency(2:end) lowerFrequency(nLinearBands)*logarithmicScale^(nLogarithmicBands+1)];
    upperFrequency = [centerFrequency(2:end) lowerFrequency(nLinearBands)*logarithmicScale^(nLogarithmicBands+2)];

    % Warp frequencies if needed
    if warped
        switch samplingFrequency
            case 8000
                    warpFactor = 0.3624;
            case 16000
                    warpFactor = 0.4595;
            case 22050
                    warpFactor = 0.52;
            case 44100
                    warpFactor = 0.6;
            otherwise
                warpFactor = 0.4595;
                warning('Using unusual sampling frequency. Do not know what to set warp-factor to.')
        end
        lowerFrequency = warpFrequency(lowerFrequency);
        centerFrequency = warpFrequency(centerFrequency);
        upperFrequency = warpFrequency(upperFrequency);
    end

    frequencies = (0:fftSize-1) * samplingFrequency/fftSize;
    melBandMatrix=zeros(nBands, fftSize);

    for iBand=1:nBands
        betweenLowerAndCenter = frequencies > lowerFrequency(iBand) & frequencies <= centerFrequency(iBand);
        betweenCenterAndUpper = frequencies < upperFrequency(iBand) & frequencies > centerFrequency(iBand);
        melBandMatrix(iBand, :) = betweenLowerAndCenter .* (frequencies - lowerFrequency(iBand)) ...
            / (centerFrequency(iBand)-lowerFrequency(iBand));
        melBandMatrix(iBand, :) = melBandMatrix(iBand, :) + betweenCenterAndUpper .* ...
            (frequencies - upperFrequency(iBand)) / (centerFrequency(iBand)-upperFrequency(iBand));
        % Normalize to have weight one
        melBandMatrix(iBand, :) = melBandMatrix(iBand, :) * 2 / (upperFrequency(iBand) - lowerFrequency(iBand));
    end

    %
    % Specify DCT matrix
    %

    % It appears that multiplication by the DCT matrix is way faster than
    % using either the dct or fft function
    dctMatrix=dctmtx(nBands);
    dctMatrix(nDctCoefficients+1:end,:)=[];

    %
    % Specify window function
    %

    if useWindow
        window=hamming(windowSize);
    else
        window=ones(windowSize,1);
    end

    
    %
    % Apply preemphasis filter if needed
    %
    
    if preemphasis
        % The Preemphasis filter is supposed to make the frequency contents of a
        % speech spectrum more white on a long term scale. (Reference needed)

        % wav = wav - 0.97*[0; wav(1:end-1)];
        % For some reason, the filter command is about twice as slow as the command above.
        % wav = filter([1 -0.97], 1, wav);

        bs=44100;
        oldtemp=0;

        for n=1:bs:length(wav)-bs
            temp = wav(n+bs-1);
            wav(n:n+bs-1) = wav(n:n+bs-1) - 0.97*[oldtemp; wav(n:n+bs-2)];
            oldtemp = temp;
        end
        if isempty(n)
            wav = wav - 0.97*[0; wav(1:end-1)];
        else
            wav(n+bs:end) = wav(n+bs:end) - 0.97*[temp; wav(n+bs:end-1)];
        end
    end

    
    %
    % Compute warped autocorrelation coefficients if needed
    %

    if warped
        warpedXcorr=zeros(mvdrOrder+1, nCols);
        circularBuffer=zeros(windowSize, mvdrOrder+1);
        oldValues=zeros(1, mvdrOrder+1);
        circularIndex=0;
        linearIndex=0;
        nMissing = windowSize - stepsize;

        for iCol=0:nCols-1
            nMissing = nMissing + stepsize;

            while nMissing > 0
                if windowSize - circularIndex > 0
                    chunksize=nMissing;
                else
                    chunksize=windowSize - circularIndex;
                end
                circularBuffer(circularIndex+1:circularIndex+chunksize,1)=wav(linearIndex+1:linearIndex+chunksize);
                linearIndex = linearIndex+chunksize;
                for lag=1:mvdrOrder
                    circularBuffer(circularIndex+1,lag+1) = warpFactor*oldValues(lag+1) - warpFactor*circularBuffer(circularIndex+1,lag) + oldValues(lag);
                    for n=circularIndex+2:circularIndex+chunksize;
                        circularBuffer(n,lag+1) = warpFactor*circularBuffer(n-1,lag+1) - warpFactor*circularBuffer(n,lag) + circularBuffer(n-1,lag);
                    end
                end
                oldValues=circularBuffer(circularIndex+chunksize, :);
                nMissing = nMissing - chunksize;
                circularIndex = mod(circularIndex + chunksize, windowSize);
            end

            for lag=0:mvdrOrder
                warpedXcorr(lag+1, iCol+1)=circularBuffer(:, 1)'*circularBuffer(:, lag+1)/windowSize;
            end
        end
    end

    %
    % Do the actual computation of the MFCCs. It gets rather ugly due to various
    % optimizations
    %

    % Constants needed in the following
    temp=find(any(melBandMatrix~=0, 1));
    first=temp(1);
    last=temp(end);
    modifiedMelBandMatrix=melBandMatrix(:,first:last);
    LOGCONST=log(10);

    if strcmp(spectrumMethod, 'vomvdr')
        mvdrOrder=round(windowSize/2);
    end

    if ~strcmp(spectrumMethod, 'fft')
        complexExponentials = exp( -2i*pi*frequencies/samplingFrequency );
        expMatrix = complexExponentials(ones(mvdrOrder+1,1), :) .^ repmat((0:mvdrOrder)', 1, fftSize);
        expMatrixNoDc = expMatrix(2:end,:);
        modifiedExpMatrix = expMatrix(:, first:last);
        modifiedExpMatrixNoDc = expMatrix(2:end,first:last);
        toeplitzIndices = toeplitz(1:mvdrOrder+1);
        toeplitzIndicesLpc = toeplitz(1:mvdrOrder);
    end

    % Pre-allocate variables
    mfcc=zeros(nDctCoefficients , nCols);
    %fprintf('Window size %d\n', size(window));

   
    switch spectrumMethod
        case 'fft'
            excerpt=zeros(fftSize,1);
            k = sum(dctMatrix(1,:)) * .5*log10(fftSize);
            for iMfcc=0:nCols-1
                excerpt(1:windowSize) = window.*wav(iMfcc*stepsize+1:iMfcc*stepsize+windowSize);
                frequencyMagnitude = abs(fft(excerpt));
                % frequencyMagnitude = abs(fft(excerpt)/sqrt(fftSize));
                if takeLog
                    melBands = log(minConst + modifiedMelBandMatrix * frequencyMagnitude(first:last))/LOGCONST;
                    mfcc(:, iMfcc+1) = dctMatrix * melBands;

                    mfcc(1, iMfcc+1) = mfcc(1, iMfcc+1) - k;
                else
                    melBands = modifiedMelBandMatrix * frequencyMagnitude(first:last);
                    mfcc(:, iMfcc+1) = dctMatrix * melBands;
                end

                %                mfcc(1, iMfcc+1) = mfcc(1, iMfcc+1) - .5*log(fftSize)/LOGCONST;
            end

        case 'orgfft'
            iMfcc = 0;
            for wavIndex=1:stepsize:length(wav)-windowSize+1
                iMfcc = iMfcc + 1;
                excerpt = window.*wav(wavIndex:wavIndex+windowSize-1);
                frequencyMagnitude = abs(fft( excerpt, fftSize));
                melBands = log10(minConst + melBandMatrix*frequencyMagnitude);
                mfcc(:, iMfcc) = dctMatrix * melBands;
            end

        case {'mvdrDirectly', 'warpedMvdrDirectly'}
            iMfcc = 0;
            for wavIndex=1:stepsize:length(wav)-windowSize+1
                iMfcc = iMfcc + 1;
                excerpt = window.*wav(wavIndex:wavIndex+windowSize-1);

                if warped
                    R=warpedXcorr(:, iMfcc);
                    Rmatrix = toeplitz(R);
                else
                    R=xcorr(excerpt, mvdrOrder, 'biased');
                    Rmatrix = toeplitz(R((end+1)/2 + (0:mvdrOrder)));
                end
                frequencyMagnitude = (1./abs(sum(conj(expMatrix).*(Rmatrix\expMatrix), 1)))';
                melBands = log10(minConst + melBandMatrix * frequencyMagnitude);
                mfcc(:, iMfcc) = dctMatrix * melBands;
            end

        case {'vomvdr','scaledVomvdr'}
            scaled=strcmp(spectrumMethod, 'scaledVomvdr')
            iMfcc = 0;
            for wavIndex=1:stepsize:length(wav)-windowSize+1
                iMfcc = iMfcc + 1;
                %                R = xcorr(excerpt, mvdrOrder, 'biased');
                %                R = R( (end+1)/2 + (0:mvdrOrder) );
                excerpt = window.*wav(wavIndex:wavIndex+windowSize-1);
                R = myxcorr(excerpt, mvdrOrder);
                
                if R(1)==0
                    mfcc(:,iMfcc)=-inf;
                    continue
                end
                
                [coef, err, reflection]  = mylevinson(R, windowSize);

                thisOrder=length(coef)-1;
%                fprintf('Variable Order MVDR order: %d\n', length(coef));

     
                mu=zeros(1,thisOrder+1);
                for k=0:thisOrder
                    for i=0:thisOrder-k
                        mu(k+1)=mu(k+1) + (thisOrder+1-k-2*i)*coef(i+1)*conj(coef(i+k+1));
                    end
                end

                reciprocalFrequencyMagnitude = mu(1) + ...
                    2*real(mu(2:end)*modifiedExpMatrixNoDc(1:thisOrder,:));
                frequencyMagnitude = sqrt((windowSize*err)./reciprocalFrequencyMagnitude);

                melBands = log(minConst + modifiedMelBandMatrix * frequencyMagnitude')/LOGCONST;
                mfcc(:, iMfcc) = dctMatrix * melBands;

            end

 
        case {'mvdr', 'warpedMvdr'}
            iMfcc = 0;
            for wavIndex=1:stepsize:length(wav)-windowSize+1
                iMfcc = iMfcc + 1;
                if warped
                    R=warpedXcorr(:, iMfcc);
                else
                    %                R = xcorr(excerpt, mvdrOrder, 'biased');
                    %                R = R( (end+1)/2 + (0:mvdrOrder) );
                    excerpt = window.*wav(wavIndex:wavIndex+windowSize-1);
                    R = myxcorr(excerpt, mvdrOrder);
                end
                if R(1)==0
                    mfcc(:,iMfcc)=-inf;
                    continue
                end
                
                [coef, err, reflection]  = levinson(R);
                if all(abs(reflection)<1) && all(isfinite(reflection)) && isfinite(err) && err > 0.001*R(1)
                    mu=zeros(1,mvdrOrder+1);
                    for k=0:mvdrOrder
                        for i=0:mvdrOrder-k
                            mu(k+1)=mu(k+1) + (mvdrOrder+1-k-2*i)*coef(i+1)*conj(coef(i+k+1));
                        end
                    end
                    % reciprocalFrequencyMagnitude = mu(1) + ...
                    %     2*real(mu(2:end)*expMatrixNoDc);
                    % frequencyMagnitude = err./reciprocalFrequencyMagnitude;
                    % melBands = log10(melBandMatrix * frequencyMagnitude');
                    % mfcc(:, iMfcc) = dctMatrix * melBands;

                    % Optimized version of the three commented lines above
                    reciprocalFrequencyMagnitude = mu(1) + ...
                        2*real(mu(2:end)*modifiedExpMatrixNoDc);
                    frequencyMagnitude = sqrt((windowSize*err)./reciprocalFrequencyMagnitude);
                else
%                    disp('Danger of numerical instability')
%                    min(find(abs(reflection)>1))
                    Rmatrix = R(toeplitzIndices);
                    frequencyMagnitude = sqrt(windowSize./abs(sum(conj(modifiedExpMatrix).*(Rmatrix\modifiedExpMatrix), 1)));
                end

                melBands = log(minConst + modifiedMelBandMatrix * frequencyMagnitude')/LOGCONST;
                mfcc(:, iMfcc) = dctMatrix * melBands;

            end

     case {'scaledMvdr'} %, 'scaledwarpedMvdr'}
            error('Sorry, but the code you wish to execute has not been finished yet.')
            iMfcc = 0;
            for wavIndex=1:stepsize:length(wav)-windowSize+1
                iMfcc = iMfcc + 1;
                if warped
                    R=warpedXcorr(:, iMfcc);
                else
                    excerpt = window.*wav(wavIndex:wavIndex+windowSize-1);
                    [R, m, idx] = myxcorr(excerpt, mvdrOrder);
                end
                if R(1)==0
                    mfcc(:,iMfcc)=-inf;
                    continue
                end
                
                [coef, err, reflection]  = levinson(R);
                if all(abs(reflection)<1) && all(isfinite(reflection)) && isfinite(err) && err > 0.001*R(1)
                    mu=zeros(1,mvdrOrder+1);
                    for k=0:mvdrOrder
                        for i=0:mvdrOrder-k
                            mu(k+1)=mu(k+1) + (mvdrOrder+1-k-2*i)*coef(i+1)*conj(coef(i+k+1));
                        end
                    end
                    % reciprocalFrequencyMagnitude = mu(1) + ...
                    %     2*real(mu(2:end)*expMatrixNoDc);
                    % frequencyMagnitude = err./reciprocalFrequencyMagnitude;

                    % Optimized version of the three commented lines above
                    reciprocalFrequencyMagnitude = mu(1) + ...
                        2*real(mu(2:end)*modifiedExpMatrixNoDc);
                    frequencyMagnitude = sqrt((windowSize*err)./reciprocalFrequencyMagnitude);
                else
%                    disp('Danger of numerical instability')
%                    min(find(abs(reflection)>1))
                    Rmatrix = R(toeplitzIndices);
                    frequencyMagnitude = sqrt(windowSize./abs(sum(conj(modifiedExpMatrix).*(Rmatrix\modifiedExpMatrix), 1)));
                end

                if warped,
                    error('Sorry, but the code you wish to execute has not been written yet.')
                end
                frequencyMagnitude = exp(log(frequencyMagnitude)*log(m)/log(frequencyMagnitude(idx-first+1)));
                
                melBands = log(minConst + modifiedMelBandMatrix * frequencyMagnitude')/LOGCONST;
                mfcc(:, iMfcc) = dctMatrix * melBands;

            end

        case {'lpc', 'warpedLpc'}
            iMfcc = 0;
            for wavIndex=1:stepsize:length(wav)-windowSize+1
                iMfcc = iMfcc + 1;
                if warped
                    R=warpedXcorr(:, iMfcc);
                else
                    % R = xcorr(excerpt, mvdrOrder, 'biased');
                    % R = R( (end+1)/2 + (0:mvdrOrder) );
                    excerpt = window.*wav(wavIndex:wavIndex+windowSize-1);
                    R = myxcorr(excerpt, mvdrOrder);
                end
                if R(1)==0
                    mfcc(:,iMfcc)=-inf;
                    continue
                end
                [coef, err, reflection] = levinson(R);
                if ~(all(abs(reflection)<1) && all(isfinite(reflection)) && isfinite(err) && err > 0.001*R(1))
%                    disp('Danger of numerical instability')
%                    min(find(abs(reflection)>1))
                    Rmatrix = R(toeplitzIndicesLpc);
                    coef = [1; -Rmatrix\R(2:end)]';
                    err = coef*toeplitz(R)*coef';
                end

                frequencyMagnitude = sqrt(err)./(abs(coef*modifiedExpMatrix))';
                melBands = log(minConst + modifiedMelBandMatrix * frequencyMagnitude)/LOGCONST;
                % frequencyMagnitude = err./(abs(coef*expMatrix).^2)';
                % melBands = log(melBandMatrix * frequencyMagnitude)/LOGCONST;
                mfcc(:, iMfcc) = dctMatrix * melBands;
            end

        otherwise
            error('Invalid spectrum method')
    end

    if ignoreSilent
        silentValue = sum(dctMatrix(1,:))*log10(minConst/sqrt(fftSize));
        silentPlaces = mfcc(1,:) < silentValue*(1-100000*eps);
        mfcc = mfcc(:, ~silentPlaces);
    end

    if nargout > 2
        % Provide debug information

        switch(spectrumMethod)
            case 'fft'
                first=1; 
                last=512;
                frequencyMagnitude=frequencyMagnitude/sqrt(fftSize);
        end

        debug.melBandMatrix=melBandMatrix;
        debug.magnitude=frequencyMagnitude;
        debug.lowerFrequency = lowerFrequency;
        debug.centerFrequency = centerFrequency;
        debug.upperFrequency = upperFrequency;
        debug.freq=(first-1:last-1)*samplingFrequency/fftSize;

        if warped
            f=[0:fftSize-1]*samplingFrequency/fftSize;
            g = warpFrequency(f);
            p=polyfit(f,g,10);
            for n=1:length(debug.freq)
                r=roots([p(1:end-1) p(end)-debug.freq(n)]);
                for m=1:length(r)
                    if isreal(r(m)) && r(m) < samplingFrequency && r(m) >= 0
                        debug.freq(n)=r(m);
                    end
                end
            end
            debug.f=f;
            debug.wf=g;
        end

        %%% DEBUG END
    end
end

