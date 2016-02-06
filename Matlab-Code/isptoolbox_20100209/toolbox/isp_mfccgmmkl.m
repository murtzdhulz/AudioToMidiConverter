%ISP_MFCCGMMKL  Definition of MFCC-GMM-Kullback-Leibler distance measure
%
% SYNTAX
%   distancemeasure = isp_tichroma
%   distancemeasure = isp_tichroma(mfccVersion)
%
% DESCRIPTION
%   Return a struct that defines a distance measure based on MFCCs, GMMs
%   and the Kullback-Leibler distance. Functions such as isp_evaluate,
%   isp_extractfeature, isp_computedistance accept this struct as input.
%
% INPUT
%   mfccVersion:
%     Either 'at' or 'vb', specifying whether to use isp_mfccat or
%     isp_mfccvb. Default: 'vb'.
%
% OUTPUT
%   distancemeasure:
%     Struct defining the distance measure.
%
% SEE ALSO
%   isp_evaluate, isp_extractfeature, isp_computedistance, isp_mfccgmmkl,
%   isp_mfccat, isp_mfccvb.
%
% HISTORY
%   Created by Jesper H. Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function output = isp_mfccgmmkl(functionName, varargin)
    if ~exist('functionName')
        functionName = 'distanceMeasureStruct';
    end

    if ischar(functionName) && ~exist(functionName)
        varargin={functionName, varargin{:}};
        functionName = 'distanceMeasureStruct';
    end

    output = feval(functionName, varargin{:});
end

function distancemeasure = distanceMeasureStruct(mfccVersion)

    if ~exist('mfccVersion')
        mfccVersion = 'at';
    end

    fs = 22050;

    % Define distance measure
    distancemeasure.name = sprintf('MFCC-GMM-KL');
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
    
    % Options for mymfcc, isp_gmmtrain and isp_computedistancematrix

    distancemeasure.options.mfccversion = mfccVersion;
    switch mfccVersion
      case 'vb'
        distancemeasure.description = sprintf('MFCC-vb+GMM+KL');
        distancemeasure.options.mfcc.samplerate = fs;
        distancemeasure.options.mfcc.nc = 20;
        distancemeasure.options.mfcc.framesize = ...
            floor(.015*distancemeasure.options.mfcc.samplerate);
        distancemeasure.options.mfcc.hopsize = ...
            floor(.0075*distancemeasure.options.mfcc.samplerate);
        [dummy, distancemeasure.options.mfcc] = ...
            isp_mfccvb([], distancemeasure.options.mfcc);
        distancemeasure.options.mfcc = ...
            rmfield(distancemeasure.options.mfcc, 'timestamp');
      case 'at'
        distancemeasure.description = sprintf('MFCC-at+GMM+KL');
        distancemeasure.options.mfcc.samplerate = fs;
        distancemeasure.options.mfcc.logAmplitude = true;
        [dummy, distancemeasure.options.mfcc] = ...
            isp_mfccat(1, distancemeasure.options.mfcc);
      case 'siggi'
        error('Not implemented yet.')
      otherwise
        error('Invalid MFCC version specified')
    end

    distancemeasure.separatechannels = false;
    [dummy, distancemeasure.options.gmm]=isp_gmmtrain([0 1], 'nMixtures', 1, ...
                                                      'diagonal', false);
    [dummy, distancemeasure.options.gmmdistance] = isp_gmmdistance(dummy, ...
                                                      dummy, ...
                                                      'silent', true, ...
                                                      'type', 'exactkl');
end

function feature = computefeature(wav, options)

    if ~iscell(wav)
        fprintf('Extracting MFCCs.\n')
        switch options.mfccversion
            case 'at'
              mfccs = isp_mfccat(wav, options.mfcc);
            case 'vb'
              mfccs = isp_mfccvb(wav, options.mfcc);
          otherwise
            error('Unknown MFCC version.')
        end
        fprintf('Training GMM.\n')
        feature = isp_gmmtrain(mfccs, options.gmm);
        fprintf('Done extracting feature.\n')
    else
        if options.gmm.nMixtures ~= length(wav)
            error('Sorry, won''t do that.')
        end
        opt = options;
        opt.gmm.nMixtures = 1;
        for i=1:length(wav)
            switch options.mfccversion
              case 'at'
                mfccs = isp_mfccat(wav{i}, options.mfcc);
              case 'vb'
                mfccs = isp_mfccvb(wav{i}, options.mfcc);
              otherwise
                error('Unknown MFCC version.')
            end
% $$$             isActive = mfccs(:,1) > -100;
% $$$             gmm(i) = isp_gmmtrain(mfccs(isActive,:)', opt);
            gmm(i) = isp_gmmtrain(mfccs, opt.gmm);
% $$$             energy(i) = mean(wav{i}.^2);
% $$$             if strcmp(options.gmmdistance.type, 'correlation')
% $$$                 covSum = 2 * gmm(i).covariance;
% $$$                 norm = 1/(sqrt(det(2*pi*covSum )));
% $$$                 gmm(i).weight = 1/sqrt(norm);
% $$$             end
        end

        feature.type = 'gmm';
        feature.weight = cat(2, gmm.weight);
% $$$         feature.weight = cat(2, gmm.weight).*energy(:)';
        feature.weight = feature.weight / sum(feature.weight);
        feature.mean = cat(2, gmm.mean);
        feature.covariance = cat(3, gmm.covariance);
        feature.nMixtures = length(wav);
        feature.diagonal = false;
    end
    
end

function distancematrix = computedistancematrix(features1, features2, options)
    distancematrix = isp_gmmdistance(features1, features2, options.gmmdistance);
end
