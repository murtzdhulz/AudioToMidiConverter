%ISP_TIRHYTHM  Define time scale insensitive measure of rhythmic distance.
%
% SYNTAX
%   distancemeasure = isp_tirhythm
%
% DESCRIPTION
%   Return a struct that defines a time scale invariant
%   distance measure based on the envelope. Functions such as isp_evaluate,
%   isp_extractfeature, isp_computedistance accept this struct as input.
%
% OUTPUT
%   distancemeasure:
%     Struct defining the distance measure. The field 'options' is yet
%     another struct that defines the behaviour of the distance measure
%     and has the following fields:
%     tMin:
%       Ignore correlation with lag smaller than this value in seconds
%       (default: 0.1).
%     tMax:
%       Ignore correlation with lag larger than this value in seconds
%       (default: 4).
%     nBands:
%       Number of bands that the autocorrelation is reduced to (default: 45)
%     variation:
%       One of the following strings specifying a variation:
%       'rawti':
%         The tempo-insensitive version using logarithmically spaced bands.
%       'linearti':
%         A non tempo-insensitive version that uses linearly spaced bands.
%       'rawseyerlehner':
%         A re-implementation of the feature from the ISMIR 2007 paper
%         "From Rhythm Patterns to Perceived Tempo" by K. Seyerlehner,
%         G. Widmer and D. Schnitzer.
%     offset:
%       (default: [ -1, 0, 1])
%     cacheRawData:
%       Boolean specifying whether the tMin, tMax, nBands and variation
%       options are applied during the feature extraction stage (if it
%       is false), thus giving a compact feature, or not until the
%       distance computation stage (if true), which makes it
%       computationally much more feasible to experiment with different
%       settings. See how to use it in isp_rhythmdemo. (default: false)
%
% EXAMPLE
%   dstMsr = isp_tirhythm;
%   %dstMsr.options.nBands = 60; % Optionally change parameters
%   file1 = fullfile(isp_toolboxpath, 'Loveshadow - The_Acorns. Seedin Time in The Oak Room - excerpt.mp3');
%   file2 = fullfile(isp_toolboxpath, 'longmidifiles', '50s Rock.mid');
%   feature1 = isp_extractfeature(file1, dstMsr);
%   feature2 = isp_extractfeature(file2, dstMsr);
%   all_features = {feature1, feature2};
%   distanceMatrix = isp_computedistance(dstMsr, all_features, all_features)
%
% SEE ALSO
%   isp_rhythmdemo, isp_evaluate, isp_extractfeature, isp_computedistance,
% isp_mfccgmmkl.
%
% HISTORY
%   2008:  Created by Jesper H. Jensen.
%   2009:  Slight restructuring and improved help text (JHJ).

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function varargout = isp_tirhythm(functionName, varargin)
    if ~exist('functionName')
        functionName = 'distanceMeasureStruct';
    end

    if ischar(functionName) && ~exist(functionName)
        varargin={functionName, varargin{:}};
        functionName = 'distanceMeasureStruct';
    end

    [varargout{1:nargout}] = feval(functionName, varargin{:});
end

function distancemeasure = distanceMeasureStruct(mfccVersion)

    fs = 8000;

    % Define distance measure
    distancemeasure.name = sprintf('Rhythm');
    distancemeasure.computefeature = ...
        ['feature = ' mfilename '(''computefeature'', ' ...
         'wav, options);'];
    distancemeasure.computedistancematrix = ...
        ['distancematrix = ' mfilename '(''computedistancematrix'', ' ...
         'features1, features2, options);'];
    distancemeasure.samplerate = fs;
    distancemeasure.mono = true;
    distancemeasure.usedFunctions = {mfilename}; % Only needed when using
                                                 % the option to compile
                                                 % with mcc.
    
    distancemeasure.options.tMin = 0.1;
    distancemeasure.options.tMax = 4;
    distancemeasure.options.nBands = 45;
    distancemeasure.options.fs = fs;
    distancemeasure.options.variation = 'rawti';
    distancemeasure.options.offset = -1:1;
    distancemeasure.options.cacheRawData = false;
end

function feature = computefeature(wav, options)
    fprintf('Extracting rhythm pattern.\n')

    [t, xcr, D] = tempo(wav, options.fs);

    feature.t = t;
    feature.xcr = xcr;

    % Copy paste from the tempo function to obtain rawxcr
    sro = 8000;
    shop = 32;
    sgsrate = sro/shop;
    acmax = round(10*sgsrate);
    mm = (mean(max(0,diff(D')')));
    eelen = length(mm);
    onsetenv = filter([1 -1], [1 -.99],mm);
    xcr2 = xcorr(onsetenv,onsetenv,acmax);
    rawxcr = xcr2(acmax+1+[0:acmax]);
    feature.rawxcr = rawxcr;

    % Version without high-pass filtering
    xcr2 = xcorr(mm,mm,acmax);
    rawxcr = xcr2(acmax+1+[0:acmax]);
    feature.rawxcrnnzm = rawxcr;

    if ~options.cacheRawData
        feature = extractData(options, feature);
    end

end


function featurevector = extractData(options, features)

    if ~iscell(features)
        features = {features};
    end

    switch options.variation
      case {'raw', 'rawti', 'rawseyerlehner'}
        tmp = cat(1, features{:});
        featurevector = cat(1, tmp.rawxcr)';
      case {'rawnz', 'rawnzti'}
        tmp = cat(1, features{:});
        featurevector = cat(1, tmp.rawxcrnnzm)';
      case {'weight', 'weightti'}
        tmp = cat(1, features{:});
        featurevector = cat(1, tmp.xcr)';

      otherwise
        error('Unknown variation.')
    end




    tMin = options.tMin;
    tMax = options.tMax;
    nBands = options.nBands;
    fs = options.fs / 32; % Frame hop size
    flen = size(featurevector, 1);
    switch options.variation
      case {'weight', 'raw', 'rawnz'}
        lags = (0:flen-1)/fs;
        validLags = (lags >= tMin) & (lags <= tMax);
        featurevector = featurevector(validLags, :);
      case 'rawseyerlehner'
        idx=[size(featurevector,1):-1:2 1:size(featurevector,1)]';
        f=[ones(10,1); zeros(length(idx)-20,1); ones(10,1)];
        
        featurevector=ifft(fft(featurevector(idx,:)).*repmat(fft(f),1,size(featurevector,2)));
        featurevector(1:(size(featurevector,1)-1)/2,:) = [];

        lags = (0:flen-1)/fs;
        validLags = (lags >= tMin) & (lags <= tMax);
        featurevector = featurevector(validLags, :);

        

      case {'weightti', 'rawti', 'rawnzti'}

        boundaries = logspace(log10(tMin), log10(tMax), nBands+2);
        % The code below is copy-paste. 'frequencies' are actually delays.
        lowerFrequency = boundaries(1:end-2);
        centerFrequency = boundaries(2:end-1);
        upperFrequency = boundaries(3:end);
        frequencies = (0:flen-1)/fs;

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

        BandMatrix = diag(1./sum(BandMatrix, 2)') * BandMatrix;

        featurevector = BandMatrix*featurevector;
      otherwise
        error('Unknown variation.')
    end

    featurevector = featurevector ./ repmat(sqrt(sum(featurevector.^2, 1)), ...
                                            size(featurevector,1), 1);
        
end


function distancematrix = computedistancematrix(features1, features2, options)

    if options.cacheRawData
        F1 = extractData(options, features1);
        F2 = extractData(options, features2);
    else
        F1 = cat(2, features1{:});
        F2 = cat(2, features2{:});
    end

    if ~isfield(options, 'offset')
        options.offset=0;
    end

    
    oF1=F1;
    oF2=F2;
    distancematrix = zeros([length(features1) length(features2) length(options.offset)]);
    for iOff=1:length(options.offset)
        %fprintf(1, 'Offsetting data by %d\n', options.offset(iOff))
        if options.offset(iOff) >= 0
            F1 = [oF1(options.offset(iOff)+1:end, :); zeros(options.offset(iOff),size(F1,2))];
        elseif options.offset(iOff) < 0
            F1 = [zeros(-options.offset(iOff),size(oF1,2)); oF1(1:end+options.offset(iOff), :); ];
            %F1 = [oF1(ones(-options.offset(iOff),1),:); F1(1:end+options.offset(iOff), :); ];
        else
            error('Invalid offset value')
        end

        distancematrix(:,:,iOff) = repmat(sum(F1.^2, 1)', 1, size(F2,2)) + ...
            repmat(sum(F2.^2, 1), size(F1,2), 1) - 2*F1'*F2;
    end

    distancematrix = min(distancematrix, [], 3);

end

function [t,xcr,D,onsetenv,sgsrate] = tempo(d,sr,tmean,tsd,onsetenv,debug)
% [t,xcr,D,onsetenv,sgsrate] = tempo(d,sr,tmean,tsd,onsetenv,debug)
%    Estimate the overall tempo of a track for the MIREX McKinney
%    contest.  
%    d is the input audio at sampling rate sr.  tmean is the mode
%    for BPM weighting (in bpm) and tsd is its spread (in octaves).
%    onsetenv is an already-calculated onset envelope (so d is
%    ignored).  debug causes a debugging plot.
%    Output t(1) is the lower BPM estimate, t(2) is the faster,
%    t(3) is the relative weight for t(1) compared to t(2).
%    xcr is the windowed autocorrelation from which the BPM peaks were picked.
%    D is the mel-freq spectrogram
%    onsetenv is the "onset strength waveform", used for beat tracking
%    sgsrate is the sampling rate of onsetenv and D.
%
% 2006-08-25 dpwe@ee.columbia.edu
% uses: localmax, fft2melmx

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

    if nargin < 3;   tmean = 120; end
    if nargin < 4;   tsd = 1.4; end
    if nargin < 5;   onsetenv = []; end
    if nargin < 6;   debug = 0; end

    sro = 8000;
    % specgram: 256 bin @ 8kHz = 32 ms / 4 ms hop
    swin = 256;
    shop = 32;
    % mel channels
    nmel = 40;
    % sample rate for specgram frames (granularity for rest of processing)
    sgsrate = sro/shop;
    % autoco out to 4 s
    acmax = round(4*sgsrate);

    D = 0;
    
    if length(onsetenv) == 0
        % no onsetenv provided - have to calculate it

        % resample to 8 kHz
        if (sr ~= sro)
            gg = gcd(sro,sr);
            d = resample(d,sro/gg,sr/gg);
            sr = sro;
        end


        %% Using spectrogram instead since specgram seems to be deprecated
        %%D = specgram(d,swin,sr,swin,swin-shop);
        %D=spectrogram(d, hanning(swin), swin-shop, swin, sr);
        %% Construct db-magnitude-mel-spectrogram
        %mlmx = fft2melmx(swin,sr,nmel);
        %D = 20*log10(max(1e-10,mlmx(:,1:(swin/2+1))*abs(D)));

        % Doing a for-loop since spectrogram use way too much memory otherwise
        D=zeros(nmel, floor(1+(length(d)-swin)/shop));
        stepsize=shop*ceil(sr*30/shop);
        mlmx = fft2melmx(swin,sr,nmel);
        idx=1;
        for n=1:stepsize:length(d)-swin+1
            tmp=spectrogram(d(n:min(length(d),n+stepsize+swin-shop)), ...
                            hanning(swin), swin-shop, swin, sr);
            % Construct db-magnitude-mel-spectrogram
            D(:,idx:idx+size(tmp,2)-1) = ...
                20*log10(max(1e-10,mlmx(:,1:(swin/2+1))*abs(tmp)));
            idx=idx+size(tmp,2);
        end            

        % Only look at the top 80 dB
        D = max(D, max(max(D))-80);

        %imgsc(D)
        
        % The raw onset decision waveform
        mm = (mean(max(0,diff(D')')));
        eelen = length(mm);

        % dc-removed mm
        onsetenv = filter([1 -1], [1 -.99],mm);

    end  % of onsetenv calc block

    % Find rough global period
    % Only use the 1st 90 sec to estimate global pd (avoid glitches?)

    maxdur = 90; % sec
    maxcol = min(round(maxdur*sgsrate),length(onsetenv));

    xcr = xcorr(onsetenv(1:maxcol),onsetenv(1:maxcol),acmax);

    % find local max in the global ac
    rawxcr = xcr(acmax+1+[0:acmax]);

    % window it around default bpm
    bpms = 60*sgsrate./([0:acmax]+0.1);
    xcrwin = exp(-.5*((log(bpms/tmean)/log(2)/tsd).^2));

    xcr = rawxcr.*xcrwin;

    xpks = localmax(xcr);  
    % will not include any peaks in first down slope (before goes below
    % zero for the first time)
    xpks(1:min(find(xcr<0))) = 0;
    % largest local max away from zero
    maxpk = max(xcr(xpks));

    % ?? then period is shortest period with a peak that approaches the max
    %maxpkthr = 0.4;
    %startpd = -1 + min(find( (xpks.*xcr) > maxpkthr*maxpk ) );
    %startpd = -1 + (find( (xpks.*xcr) > maxpkthr*maxpk ) );

    % no, just largest peak after windowing
    startpd = -1 + find((xpks.*xcr) == max(xpks.*xcr));

    % ??Choose acceptable peak closest to 120 bpm
    %[vv,spix] = min(abs(60./(startpd/sgsrate) - 120));
    %startpd = startpd(spix);
    % No, just choose shortest acceptable peak
    startpd = startpd(1);

    t = 60/(startpd/sgsrate);

    % Choose best peak out of .33 .5 2 3 x this period
    candpds = round([.33 .5 2 3]*startpd);
    candpds = candpds(candpds < acmax);

    [vv,xx] = max(xcr(1+candpds));

    startpd2 = candpds(xx);
    vvm = xcr(1+startpd);
    pratio = vvm/(vvm+vv);

    t = [60/(startpd/sgsrate) 60/(startpd2/sgsrate) pratio];

    % ensure results are lowest-first
    if t(2) < t(1)
        t([1 2]) = t([2 1]);
        t(3) = 1-t(3);
    end  

    startpd = (60/t(1))*sgsrate;
    startpd2 = (60/t(2))*sgsrate;

    %  figure
    %  disp(['tmean=',num2str(tmean),' tsd=',num2str(tsd),' maxpk=',num2str(startpd)]);
    %  subplot(211)
    %  plot([0:acmax],xcrwin/max(abs(xcrwin)),[0:acmax],xcr/max(abs(xcr)),...
    %       [startpd startpd],[-1 1],'-r',[startpd2 startpd2],[-1 1],'-c')
    %  subplot(212)
    %  bpms(1) = bpms(2);
    %  plot(bpms,xcrwin/max(abs(xcrwin)),bpms,xcr/max(abs(xcr)),...
    %       [t(1) t(1)],[-1 1],'-r',[t(2) t(2)],[-1 1],'-c')

    if debug > 0

        % Report results and plot weighted autocorrelation with picked peaks
        disp(['Global bt pd = ',num2str(t(1)),' @ ',num2str(t(3)),' / ',num2str(t(2)),' bpm']);

        subplot(414)
        plot([0:acmax],xcr,'-b', ...
             [0:acmax],xcrwin*maxpk,'-r', ...
             [startpd startpd], [min(xcr) max(xcr)], '-g', ...
             [startpd2 startpd2], [min(xcr) max(xcr)], '-c');
        grid;

    end
end
% Read in all the tempo settings
% for i = 1:20; f = fopen(['mirex-beattrack/train/train',num2str(i),'-tempo.txt']); r(i,:) = fscanf(f, '%f\n'); fclose(f); end
function m = localmax(x)
% return 1 where there are local maxima in x (columnwise).
% don't include first point, maybe last point

    [nr,nc] = size(x);

    if nr == 1
        lx = nc;
    elseif nc == 1
        lx = nr;
        x = x';
    else
        lx = nr;
    end

    if (nr == 1) || (nc == 1)

        m = (x > [x(1),x(1:(lx-1))]) & (x >= [x(2:lx),1+x(lx)]);

        if nc == 1
            % retranspose
            m = m';
        end
        
    else
        % matrix
        lx = nr;
        m = (x > [x(1,:);x(1:(lx-1),:)]) & (x >= [x(2:lx,:);1+x(lx,:)]);
    end
end

function [wts,binfrqs] = fft2melmx(nfft, sr, nfilts, width, minfrq, maxfrq, htkmel, constamp)
% wts = fft2melmx(nfft, sr, nfilts, width, minfrq, maxfrq, htkmel, constamp)
%      Generate a matrix of weights to combine FFT bins into Mel
%      bins.  nfft defines the source FFT size at sampling rate sr.
%      Optional nfilts specifies the number of output bands required 
%      (else one per bark), and width is the constant width of each 
%      band relative to standard Mel (default 1).
%      While wts has nfft columns, the second half are all zero. 
%      Hence, Mel spectrum is fft2melmx(nfft,sr)*abs(fft(xincols,nfft));
%      minfrq is the frequency (in Hz) of the lowest band edge;
%      default is 0, but 133.33 is a common standard (to skip LF).
%      maxfrq is frequency in Hz of upper edge; default sr/2.
%      You can exactly duplicate the mel matrix in Slaney's mfcc.m
%      as fft2melmx(512, 8000, 40, 1, 133.33, 6855.5, 0);
%      htkmel=1 means use HTK's version of the mel curve, not Slaney's.
%      constamp=1 means make integration windows peak at 1, not sum to 1.
% 2004-09-05  dpwe@ee.columbia.edu  based on fft2barkmx

    if nargin < 2;     sr = 8000;      end
    if nargin < 3;     nfilts = 40;    end
    if nargin < 4;     width = 1.0;    end
    if nargin < 5;     minfrq = 0;     end  % default bottom edge at 0
    if nargin < 6;     maxfrq = sr/2;  end  % default top edge at nyquist
    if nargin < 7;     htkmel = 0;     end
    if nargin < 8;     constamp = 0;   end


    wts = zeros(nfilts, nfft);

    % Center freqs of each FFT bin
    fftfrqs = [0:(nfft/2)]/nfft*sr;

    % 'Center freqs' of mel bands - uniformly spaced between limits
    minmel = hz2mel(minfrq, htkmel);
    maxmel = hz2mel(maxfrq, htkmel);
    binfrqs = mel2hz(minmel+[0:(nfilts+1)]/(nfilts+1)*(maxmel-minmel), htkmel);

    binbin = round(binfrqs/sr*(nfft-1));

    for i = 1:nfilts
        %  fs = mel2hz(i + [-1 0 1], htkmel);
        fs = binfrqs(i+[0 1 2]);
        % scale by width
        fs = fs(2)+width*(fs - fs(2));
        % lower and upper slopes for all bins
        loslope = (fftfrqs - fs(1))/(fs(2) - fs(1));
        hislope = (fs(3) - fftfrqs)/(fs(3) - fs(2));
        % .. then intersect them with each other and zero
        %  wts(i,:) = 2/(fs(3)-fs(1))*max(0,min(loslope, hislope));
        wts(i,1+[0:(nfft/2)]) = max(0,min(loslope, hislope));

        % actual algo and weighting in feacalc (more or less)
        %  wts(i,:) = 0;
        %  ww = binbin(i+2)-binbin(i);
        %  usl = binbin(i+1)-binbin(i);
        %  wts(i,1+binbin(i)+[1:usl]) = 2/ww * [1:usl]/usl;
        %  dsl = binbin(i+2)-binbin(i+1);
        %  wts(i,1+binbin(i+1)+[1:(dsl-1)]) = 2/ww * [(dsl-1):-1:1]/dsl;
        % need to disable weighting below if you use this one

    end

    if (constamp == 0)
        % Slaney-style mel is scaled to be approx constant E per channel
        wts = diag(2./(binfrqs(2+[1:nfilts])-binfrqs(1:nfilts)))*wts;
    end

    % Make sure 2nd half of FFT is zero
    wts(:,(nfft/2+1):nfft) = 0;
    % seems like a good idea to avoid aliasing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = mel2hz(z, htk)
%   f = mel2hz(z, htk)
%   Convert 'mel scale' frequencies into Hz
%   Optional htk = 1 means use the HTK formula
%   else use the formula from Slaney's mfcc.m
% 2005-04-19 dpwe@ee.columbia.edu

    if nargin < 2
        htk = 0;
    end

    if htk == 1
        f = 700*(10.^(z/2595)-1);
    else
        
        f_0 = 0; % 133.33333;
        f_sp = 200/3; % 66.66667;
        brkfrq = 1000;
        brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
        logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)

        linpts = (z < brkpt);

        f = 0*z;

        % fill in parts separately
        f(linpts) = f_0 + f_sp*z(linpts);
        f(~linpts) = brkfrq*exp(log(logstep)*(z(~linpts)-brkpt));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = hz2mel(f,htk)
%  z = hz2mel(f,htk)
%  Convert frequencies f (in Hz) to mel 'scale'.
%  Optional htk = 1 uses the mel axis defined in the HTKBook
%  otherwise use Slaney's formula
% 2005-04-19 dpwe@ee.columbia.edu

    if nargin < 2
        htk = 0;
    end

    if htk == 1
        z = 2595 * log10(1+f/700);
    else
        % Mel fn to match Slaney's Auditory Toolbox mfcc.m

        f_0 = 0; % 133.33333;
        f_sp = 200/3; % 66.66667;
        brkfrq = 1000;
        brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
        logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)

        linpts = (f < brkfrq);

        z = 0*f;

        % fill in parts separately
        z(linpts) = (f(linpts) - f_0)/f_sp;
        z(~linpts) = brkpt+(log(f(~linpts)/brkfrq))./log(logstep);

    end
end