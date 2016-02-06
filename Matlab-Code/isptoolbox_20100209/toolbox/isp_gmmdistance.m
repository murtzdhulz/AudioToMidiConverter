%ISP_GMMDISTANCE  Compute distance between two Gaussian mixture models.
%
% SYNTAX
%   [distancematrix, options] = isp_computedistancematrix(features1, features2, options ...)
%
% DESCRIPTION
%   Computes the distance between Gaussian mixture models.
%
% INPUT
%   features1, features2:
%     Cell array of gmmstructs as returned by isp_gmmtrain.
%   options ...:
%     The following parameters can be set as field names in structs, or
%     be specified as field/value pairs:
%     type:
%       A string indicating how to compute the distance. Possibilities:
%       'montecarlokl':
%         Kullback-Leibler distance computed by Monte Carlo sampling.
%       'fastkl':
%         Calculate a fast, but less accurate approximation of the KL-distance.
%       'exactkl':
%         Exact KL-distance. Only possible with a single mixture.
%       'euclidean':
%         The L2 distance between the GMMs.
%       'correlation':
%         Inner product between the GMM's probability density
%         functions normalized to unit L2 norm, i.e., the correlation
%         coefficient or a continuous version of the cosine between
%         vectors. The actual output is one minut the inner product to
%         ensure that smaller is better.
%       A prefix of 'emd', e.g. 'emdexact', computes earth movers
%       distance (EMD) using the specified distance as cost. You need
%       Simon Dixons Matlab wrapper for Yossi Rubners EMD implementation.
%       Default is 'exactkl' for a single Gaussian and 'fastkl' for mixtures.
%     nSamples:
%       An integer specifying the number of samples when 'montecarlokl' is
%       chosen as type.
%       Default is 200.
%     weightThreshold:
%       When doing Monte Carlo sampling, do not generate random samples from
%       Gaussian mixtures with weight less than this value.
%       Default is 0.01.
%     silent:
%       Boolean specifying if information about the progress shall be
%       suppresed or not. Default: false.
%     forceDiagonal:
%       Boolean specifying if Gaussians with full covariance matrices shall be
%       converted to diagonal matrices. Provided for test purposes.
%       true: Convert full covariance matrices to diagonal.
%       false: Leave covariance matrices untouched. 
%       Default is false.
%     cache1, cache2:
%       Caching of e.g. matrix inverses for 'features1' and 'features2',
%       respectively. In reality, isp_gmmdistance not only return two
%       parameters, but actually return [distancematrix, options, cache1,
%       cache2]. If e.g. features2 never change, one can speed things up
%       by giving the returned cache2 as input. This feture is not
%       heavily tested.
%
% OUTPUT
%   Matrix of distances. Element (i,j) is the distance between the GMMs in
%   features1{i} and features2{j}.
%
% EXAMPLE
%     distancematrix = isp_computedistancematrix(features)
%     distancematrix = isp_computedistancematrix({gmm1, gmm2, gmm3}, options)
%
% SEE ALSO
%   isp_gmmprobability, isp_gmmrand, isp_gmmtrain.
%
% HISTORY
%   2006-03-31: JHJ - First version for the ISP toolbox.
%   2006-05-07: JHJ - Modified to accept two cell arrays of feature vectors
%   2007-01-22: JHJ - Added 'correlation' distance type
%   2007-04-01: JHJ - Added 'emd' prefix and added optimized code for
%                    'correlation' for diagonal covariance matrices.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.



function [distancematrix, options, cache1, cache2]=isp_gmmdistance(features1, features2, varargin)

    if nargin < 2
        error('Not enough input arguments.')
    end
    
    if ~iscell(features1)
        features1={features1};
    end

    if ~iscell(features2)
        features2={features2};
    end

    if ~(strcmp(features1{1}.type, 'gmm') && strcmp(features2{1}.type, 'gmm'))
        error('Inputs are not Gaussian mixture models.')
    end

    twoSets = ~isequal(features1, features2);

    options=struct('nSamples', 200, ...
                   'weightThreshold', 0.01, ...
                   'silent', false, ...
                   'type', '', ...
                   'forceDiagonal', false, ...
                   'euclideanweight', [], ...
                   'cache1', struct, ...
                   'cache2', struct);

    options = isp_interpretarguments(options, varargin{:});
    
    nSongs1=length(features1);
    nSongs2=length(features2);
    nDim=size(features1{1}.mean,1);
    nMixtures1=features1{1}.nMixtures;
    nMixtures2=features2{1}.nMixtures;
    mcSamples = options.nSamples;

    cache1 = options.cache1;
    cache2 = options.cache2;

    if isempty(options.type)
        if nMixtures1 > 1 || nMixtures2 > 1
            options.type = 'fastkl';
        else
            options.type = 'exactkl';
        end
    end

    if options.forceDiagonal
        if ~features1{1}.diagonal
            for n=1:nSongs1
                features1{n}.covariance = diag(features1{n}.covariance);
                features1{n}.diagonal = true;
            end
        end
        
        if ~features2{1}.diagonal
            for n=1:nSongs2
                features2{n}.covariance = diag(features2{n}.covariance);
                features2{n}.diagonal = true;
            end
        end
    end
    
    % This ought to be done in the non-diagonal case too
    if features1{1}.diagonal
        for n=1:nSongs1
            zeroVar = any( features1{n}.covariance<=0, 1);
            if any(zeroVar)
                features1{n}.weight(:,zeroVar) = 0;
                features1{n}.covariance(:,zeroVar) = 1;
            end
        end
    end

    % This ought to be done in the non-diagonal case too
    if features2{1}.diagonal
        for n=1:nSongs2
            zeroVar = any( features2{n}.covariance<=0, 1);
            if any(zeroVar)
                features2{n}.weight(:,zeroVar) = 0;
                features2{n}.covariance(:,zeroVar) = 1;
            end
        end
    end

    type=options.type;
    wThreshold = options.weightThreshold;
    silent = options.silent;

    
    
    if ~isempty(regexp(type, '^emd'))
        if ~exist('emd_wrapper', 'file')
            error(['Cannot find ''emd_wrapper''. You must obtain and compile ' ...
                   'this file separately.'])
        end
        
        for iSong1=1:nSongs1
            for iMixture1=1:nMixtures1
                single1{iSong1}{iMixture1}.nMixtures = 1;
                single1{iSong1}{iMixture1}.type = 'gmm';
                single1{iSong1}{iMixture1}.weight = 1;
                single1{iSong1}{iMixture1}.diagonal = features1{iSong1}.diagonal;
                single1{iSong1}{iMixture1}.mean = features1{iSong1}.mean(:,iMixture1);
                if single1{iSong1}{iMixture1}.diagonal
                    single1{iSong1}{iMixture1}.covariance = features1{iSong1}.covariance(:,iMixture1);
                else
                    single1{iSong1}{iMixture1}.covariance = features1{iSong1}.covariance(:,:,iMixture1);
                end
                weights1{iSong1} = features1{iSong1}.weight;
            end
        end        

        for iSong2=1:nSongs2
            for iMixture2=1:nMixtures2
                single2{iSong2}{iMixture2}.nMixtures = 1;
                single2{iSong2}{iMixture2}.type = 'gmm';
                single2{iSong2}{iMixture2}.weight = 1;
                single2{iSong2}{iMixture2}.diagonal = features2{iSong2}.diagonal;
                single2{iSong2}{iMixture2}.mean = features2{iSong2}.mean(:,iMixture2);
                if single2{iSong2}{iMixture2}.diagonal
                    single2{iSong2}{iMixture2}.covariance = features2{iSong2}.covariance(:,iMixture2);
                else
                    single2{iSong2}{iMixture2}.covariance = features2{iSong2}.covariance(:,:,iMixture2);
                end
                weights2{iSong2} = features2{iSong2}.weight;
            end
        end        

        emdoptions = options;
        emdoptions.type = ...
            regexprep(emdoptions.type, '^emd', '')

        distancematrix=zeros(nSongs1, nSongs2);
        song1ChunkSize = 30;
        for iSong1=1:song1ChunkSize:nSongs1
            iSong1
            song1Idx = iSong1:min(nSongs1, iSong1+song1ChunkSize-1);
            dstMtx = isp_computedistancematrix(cat(2, single1{song1Idx}), ...
                                               cat(2, single2{:}), ...
                                               emdoptions);
            
            idx1=1;
            for iSong2 = 1:length(song1Idx)
                nSingle1Mixtures = length(single1{song1Idx(iSong2)});
                idx2=1;
                for jSong2=1:nSongs2
                    nSingle2Mixtures = length(single2{jSong2});
                    distancematrix(song1Idx(iSong2), jSong2) = ...
                        emd_wrapper(dstMtx(idx1:idx1+nSingle1Mixtures-1,idx2:idx2+nSingle2Mixtures-1), weights1{iSong1}, weights2{jSong2});
                    idx2 = idx2 + nSingle2Mixtures;
                end
                idx1 = idx1 + nSingle1Mixtures;
            end
        end
        return

    end




    if ~isempty(regexp(type, '^weighted'))
        %fprintf('Weighting features')
        type = regexprep(type, '^weighted', '');
        W = options.euclideanweight;

        if (features1{n}.diagonal || features2{n}.diagonal) && ~isvector(W)
            error(['Sorry, linearly transforming diagonal covariance ' ...
                   'matrices by a non-diagonal matrix is not implemented.'])
        end

        if isvector(W)
            W = diag(W);
        end

        for n=1:nSongs1
            features1{n}.mean = W * features1{n}.mean;
            if features1{n}.diagonal
                features1{n}.covariance = W * W * features1{n}.covariance;
            else            
                for m=1:size(features1{n}.covariance, 3)
                    features1{n}.covariance(:,:,m) = ...
                        W * features1{n}.covariance(:,:,m) * W';
                end
            end
        end
        
        for n=1:nSongs2
            features2{n}.mean = W * features2{n}.mean;
            if features2{n}.diagonal
                features2{n}.covariance = W * W * features2{n}.covariance;
            else            
                for m=1:size(features2{n}.covariance, 3)
                    features2{n}.covariance(:,:,m) = ...
                        W * features2{n}.covariance(:,:,m) * W';
                end
            end
        end

    end
    
    
    
    
    
    switch type
      case 'montecarlokl'
        distancematrix=zeros(nSongs2, nSongs1);
        distancematrix2=zeros(nSongs1, nSongs2);
        mfccVector=zeros(nDim, mcSamples*nSongs1);
        probs=zeros(mcSamples*nSongs1,1);
        for iSong=1:nSongs1
            g=features1{iSong};
            g.weight(g.weight<wThreshold )=0;
            mfccVector(:, (iSong-1)*mcSamples+1:iSong*mcSamples)=isp_gmmrand(g, mcSamples, true);
        end

        if twoSets
            % Calculate distance to each song
            selfProb1 = zeros(nSongs1, 1);
            selfProb2 = zeros(nSongs2, 1);

            for iSong=1:nSongs1
                selfProb1(iSong) = mean(log( isp_gmmprobability(features1{iSong}, ...
                                                                mfccVector(:, (iSong-1)*mcSamples+1:iSong*mcSamples) ) ));
            end

            for iSong=1:nSongs2
                if ~silent, fprintf('Song %d of %d\n', iSong, nSongs2); end
                probs=log(isp_gmmprobability(features2{iSong}, mfccVector));
                distancematrix2(:,iSong)=mean(reshape(probs, mcSamples, nSongs1),1);
            end

            mfccVector=zeros(nDim, mcSamples*nSongs2);
            probs=zeros(mcSamples*nSongs2,1);
            for iSong=1:nSongs2
                g=features2{iSong};
                g.weight(g.weight<wThreshold )=0;
                mfccVector(:, (iSong-1)*mcSamples+1:iSong*mcSamples)=isp_gmmrand(g, mcSamples, true);
                selfProb2(iSong) = mean(log( isp_gmmprobability(features2{iSong}, ...
                                                                mfccVector(:, (iSong-1)*mcSamples+1:iSong*mcSamples) ) ));
            end
        end

        % Calculate distance to each song
        for iSong=1:nSongs1
            if ~silent, fprintf('Song %d of %d\n', iSong, nSongs1); end
            probs=log(isp_gmmprobability(features1{iSong}, mfccVector));
            distancematrix(:,iSong)=mean(reshape(probs, mcSamples, nSongs2),1);
        end

        % Normalize
        if twoSets
            distancematrix = - distancematrix' - distancematrix2 ...
                + repmat(selfProb1,1,nSongs2) ...
                + repmat(selfProb2',nSongs1,1);
        else
            distancematrix = - distancematrix-distancematrix' ...
                + repmat(diag(distancematrix),1,nSongs1) ...
                + repmat(diag(distancematrix)',nSongs1,1);
        end
    

        
      case 'fastkl'
        distancematrix=zeros(nSongs2, nSongs1);
        distancematrix2=zeros(nSongs1, nSongs2);
        nSamples=nMixtures1;
        mfccVector=zeros(nDim, nSamples*nSongs1);
        probs=zeros(nSamples*nSongs1,1);
        weightVector=zeros(size(probs));
        for iSong=1:nSongs1
            g=features1{iSong};
            m=g.mean;
            w=g.weight;
            mfccVector(:, (iSong-1)*nSamples+1:iSong*nSamples)=m;
            weightVector((iSong-1)*nSamples+1:iSong*nSamples)=w;
        end
        weightMatrix=reshape(weightVector, nSamples, nSongs1);
        
        % Calculate distance to each song
        if twoSets
            selfProb1 = zeros(nSongs1, 1);
            selfProb2 = zeros(nSongs2, 1);

            for iSong=1:nSongs1
                selfProb1(iSong) = features1{iSong}.weight*log( isp_gmmprobability(features1{iSong}, ...
                                                                  mfccVector(:, (iSong-1)*nSamples+1:iSong*nSamples) ) );
            end

            for iSong=1:nSongs2
                if ~silent, fprintf('Song %d of %d\n', iSong, nSongs2); end
                probs=log(isp_gmmprobability(features2{iSong}, mfccVector));
                distancematrix2(:,iSong)=sum(reshape(probs, nSamples, nSongs1).*weightMatrix,1)';
            end

            nSamples=nMixtures2;
            mfccVector=zeros(nDim, nSamples*nSongs2);
            probs=zeros(nSamples*nSongs2,1);
            weightVector=zeros(size(probs));

            for iSong=1:nSongs2
                g=features2{iSong};
                m=g.mean;
                w=g.weight;
                mfccVector(:, (iSong-1)*nSamples+1:iSong*nSamples)=m;
                weightVector((iSong-1)*nSamples+1:iSong*nSamples)=w;
                selfProb2(iSong) = w*log( isp_gmmprobability(features2{iSong}, ...
                                                             mfccVector(:, (iSong-1)*nSamples+1:iSong*nSamples) ) );
                
            end
            weightMatrix=reshape(weightVector, nSamples, nSongs2);
        end
        % Calculate distance to each song
        for iSong=1:nSongs1
            if ~silent, fprintf('Song %d of %d\n', iSong, nSongs1); end
            probs=log(isp_gmmprobability(features1{iSong}, mfccVector));
            distancematrix(:,iSong)=sum(reshape(probs, nSamples, nSongs2).*weightMatrix,1)';
        end
        % Normalize
        if twoSets
            distancematrix = - distancematrix' - distancematrix2 ...
                + repmat(selfProb1,1,nSongs2) ...
                + repmat(selfProb2',nSongs1,1);
        else
            distancematrix = - distancematrix-distancematrix' ...
                + repmat(diag(distancematrix),1,nSongs1) ...
                + repmat(diag(distancematrix)',nSongs1,1);
        end
    

      case {'correlation', 'euclidean'}
        distancematrix=zeros(nSongs1, nSongs2);

        euclidean = strcmp(type, 'euclidean');

        if ~twoSets
            features2 = features1;
            nMixtures2 = nMixtures1;
        end
        
        if euclidean
            norm1 = zeros(nSongs1, 1);
            norm2 = zeros(nSongs2, 1);
        end

        if features1{1}.diagonal && features2{1}.diagonal
            % Diagonal covariance matrices with 1 mixture

            % Normalize to unit norm
            for iSong1 = 1:nSongs1
                weight = features1{iSong1}.weight;
                cov = features1{iSong1}.covariance;
                means = features1{iSong1}.mean;
                norm = 0;
                for kMix1 = 1:nMixtures1
                    for lMix1 = 1:nMixtures1
                        covSum = cov(:,kMix1) + cov(:,lMix1);
                        meanDiff = (means(:,kMix1)-means(:,lMix1));
                        norm = norm ...
                               +  weight(kMix1) * weight(lMix1) ...
                               * 1/(sqrt(prod(2*pi*covSum ))) ...
                               * exp( -.5*( (meanDiff./covSum).'*meanDiff) );
                    end
                end
                if euclidean
                    norm1(iSong1) = norm;
                else
                    features1{iSong1}.weight = weight/sqrt(norm);
                end
            end
            
            % Normalize to unit norm
            for iSong2 = 1:nSongs2
                weight = features2{iSong2}.weight;
                cov = features2{iSong2}.covariance;
                means = features2{iSong2}.mean;
                norm = 0;
                for kMix2 = 1:nMixtures2
                    for lMix2 = 1:nMixtures2
                        covSum = (cov(:,kMix2) + cov(:,lMix2));
                        meanDiff = (means(:,kMix2)-means(:,lMix2));
                        norm = norm ...
                               +  weight(kMix2) * weight(lMix2) ...
                               * 1/(sqrt(prod(2*pi*covSum ))) ...
                               * exp( -.5*( (meanDiff./covSum).'*meanDiff) );
                    end
                end
                if euclidean
                    norm2(iSong2) = norm;
                else
                    features2{iSong2}.weight = weight/sqrt(norm);
                end
            end

            if nMixtures1==1 && nMixtures2==1
                temp = cat(1, features1{:});
                cov1 = cat(2, temp.covariance);
                mean1 = cat(2, temp.mean);
                weight1 = cat(2, temp.weight);
                
                temp = cat(1, features2{:});
                cov2 = cat(2, temp.covariance);
                mean2 = cat(2, temp.mean);
                weight2 = cat(2, temp.weight);

                for iSong1 = 1:nSongs1
                    covSum = cov1(:,iSong1*ones(nSongs2,1)) + cov2;
                    meanDiff = mean1(:,iSong1*ones(nSongs2,1)) - mean2;
                    distancematrix(iSong1, :) = ...
                        weight1(iSong1) * (weight2 ...
                                           ./(sqrt(prod(2*pi*covSum,1))) ...
                                           .* exp( -.5*sum((meanDiff.^2)./covSum,1)));
                end

            else


                mix1indices = repmat(1:nMixtures1, 1, nMixtures2);
                mix2indices = reshape(repmat(1:nMixtures1, nMixtures2, 1), 1, ...
                                      nMixtures1*nMixtures2);
                for iSong1 = 1:nSongs1
                    weight1 = features1{iSong1}.weight;
                    cov1 = features1{iSong1}.covariance;
                    mean1 = features1{iSong1}.mean;
                    for jSong2 = 1:nSongs2
                        weight2 = features2{jSong2}.weight;
                        cov2 = features2{jSong2}.covariance;
                        mean2 = features2{jSong2}.mean;

                        covSum = cov1(:, mix1indices) + cov2(:, mix2indices);
                        meanDiff = mean1(:, mix1indices) - mean2(:, mix2indices);
                        distancematrix(iSong1, jSong2) = sum( ...
                            weight1(mix1indices).* weight2(mix2indices) ...
                            ./sqrt(prod(2*pi*covSum, 1)) ...
                            .* exp( -.5*sum( meanDiff.^2./covSum, 1)) ...
                            );
                    end
                end
            end

        else

            % Make sure both feature sets have full covariance matrices
            if features1{1}.diagonal
                covariance = zeros([nDim nDim nMixtures1]);
                for n=1:nSongs1
                    for m=1:nMixtures1
                        covariance(:,:,m) = diag(features1{n}.covariance(:, m));
                    end
                    features1{n}.covariance = covariance;
                end
            end
            if features2{1}.diagonal
                covariance = zeros([nDim nDim nMixtures2]);
                for n=1:nSongs2
                    for m=1:nMixtures2
                        covariance(:,:,m) = diag(features2{n}.covariance(:, m));
                    end
                    features2{n}.covariance = covariance;
                end
            end

            % Compute norm (correlation: also normalize to unit norm)
            for iSong1 = 1:nSongs1
                weight = features1{iSong1}.weight;
                cov = features1{iSong1}.covariance;
                means = features1{iSong1}.mean;
                norm = 0;
                for kMix1 = 1:nMixtures1
                    for lMix1 = 1:nMixtures1
                        covSum = cov(:,:,kMix1) + cov(:,:,lMix1);
                        meanDiff = (means(:,kMix1)-means(:,lMix1));
                        norm = norm ...
                               +  weight(kMix1) * weight(lMix1) ...
                               * 1/(sqrt(det(2*pi*covSum ))) ...
                               * exp( -.5*( meanDiff.'*(covSum\meanDiff)) );
                    end
                end
                if euclidean
                    norm1(iSong1) = norm;
                else
                    features1{iSong1}.weight = weight/sqrt(norm);
                end
            end
            
            % Compute norm (correlation: also normalize to unit norm)
            for iSong2 = 1:nSongs2
                weight = features2{iSong2}.weight;
                cov = features2{iSong2}.covariance;
                means = features2{iSong2}.mean;
                norm = 0;
                for kMix2 = 1:nMixtures2
                    for lMix2 = 1:nMixtures2
                        covSum = (cov(:,:,kMix2) + cov(:,:,lMix2));
                        meanDiff = (means(:,kMix2)-means(:,lMix2));
                        norm = norm ...
                               +  weight(kMix2) * weight(lMix2) ...
                               * 1/(sqrt(det(2*pi*covSum ))) ...
                               * exp( -.5*( meanDiff.'*(covSum\meanDiff)) );
                    end
                end
                if euclidean
                    norm2(iSong2) = norm;
                else
                    features2{iSong2}.weight = weight/sqrt(norm);
                end
            end

            for iSong1 = 1:nSongs1
                weight1 = features1{iSong1}.weight;
                cov1 = features1{iSong1}.covariance;
                mean1 = features1{iSong1}.mean;
                for jSong2 = 1:nSongs2
                    weight2 = features2{jSong2}.weight;
                    cov2 = features2{jSong2}.covariance;
                    mean2 = features2{jSong2}.mean;
                    for kMix1 = 1:nMixtures1
                        for lMix2 = 1:nMixtures2
                            covSum = cov1(:,:,kMix1) + cov2(:,:,lMix2);
                            meanDiff = (mean1(:,kMix1)-mean2(:,lMix2));
                            distancematrix(iSong1, jSong2) = ...
                                distancematrix(iSong1, jSong2) ...
                                +  weight1(kMix1) * weight2(lMix2) ...
                                * 1/(sqrt(det(2*pi*covSum ))) ...
                                * exp( -.5*( meanDiff.'*(covSum\meanDiff)) );
                        end
                    end
                end
            end
        end

        if euclidean
            % Add the squared norms to the double of the product
            distancematrix = repmat(norm1(:), 1, nSongs2) ...
                + repmat(norm2(:)', nSongs1, 1) - 2*distancematrix;
        else
            % Create distance measure where 0 is equality and 1 is uncorrelated
            distancematrix = 1 - distancematrix;
        end
        
        

        
      case 'exactkl'
        if nMixtures1~=1 || ( twoSets && nMixtures2~=1 )
            error(['Unable to compute exact KL-distance for other ' ...
                   'numbers of mixtures than 1.'])
        end

        if ~twoSets
            features2 = features1;
        end

        % Things could be done much more efficiently for diagonal covariance
        % matrices, but for now, they are converted to full matrices.

        if features1{1}.diagonal && features2{1}.diagonal

            % Precompute stuff
            mean1 = zeros([nDim nSongs1]);
            
            for iSong=1:nSongs1
                mean1(:,iSong) = features1{iSong}.mean;
            end

            % Compute non-symmetric KL-distance if two sets of feature
            % vectors were passed as arguments
            if twoSets
                mean2 = zeros([nDim nSongs2]);
                for iSong=1:nSongs2
                    mean2(:,iSong) = features2{iSong}.mean;
                end

                invCov1 = zeros([nDim nSongs1]);
                for iSong=1:nSongs1
                    invCov1(:,iSong) = 1./features1{iSong}.covariance;
                end

                distancematrix2 = zeros(nSongs1, nSongs2);
                intermediateResult1 = sum(mean1.^2.*invCov1, 1).';
                intermediateResult2 = mean1.*invCov1;
                
                for iSong=1:nSongs2
                    igmmCov = features2{iSong}.covariance;
                    distancematrix2(:, iSong) = ...
                        invCov1.'*igmmCov + ...
                        invCov1.'*(mean2(:, iSong).^2) + ...
                        intermediateResult1 - ...
                        intermediateResult2.'*(2*mean2(:, iSong)) - ...
                        nDim;
                end
                distancematrix2 = .5*distancematrix2;
            else
                mean2=mean1;
            end

            % Compute non-symmetric KL-distance
            invCov2 = zeros([nDim nSongs2]);
            for iSong=1:nSongs2
                invCov2(:,iSong) = 1./features2{iSong}.covariance;
            end
            
            distancematrix = zeros(nSongs1, nSongs2);
            intermediateResult1 = sum(mean2.^2.*invCov2, 1).';
            intermediateResult2 = mean2.*invCov2;

            for iSong=1:nSongs1
                if ~silent, fprintf('Song %d of %d\n', iSong, nSongs1); end
                igmmCov = features1{iSong}.covariance;

                distancematrix(:, iSong) = ...
                    invCov2.'*igmmCov + ...
                    invCov2.'*(mean1(:, iSong).^2) + ...
                    intermediateResult1 - ...
                    intermediateResult2.'*(2*mean1(:, iSong)) - ...
                    nDim;
                
            end
            distancematrix = .5*distancematrix;
        else
            % Full covariance matrix

            % Handle the very special case when only one of the feature sets have full
            % covariance matrix.
            if features1{1}.diagonal
                for n=1:nSongs1
                    features1{n}.covariance = diag(features1{n}.covariance);
                end
            end

            if features2{1}.diagonal
                for n=1:nSongs2
                    features2{n}.covariance = diag(features2{n}.covariance);
                end
            end

            % Precompute stuff
            if isfield(cache1, 'mean1')
                mean1 = cache1.mean1;
            else
                mean1 = zeros([nDim nSongs1]);
                for iSong1=1:nSongs1
                    mean1(:,iSong1) = features1{iSong1}.mean;
                end
            end

            if isfield(cache2, 'invCov2')
                invCov2 = cache2.invCov2;
            else
                invCov2 = zeros([nDim nDim nSongs2]);
                for iSong2=1:nSongs2
                    invCov2(:,:,iSong2) = inv(features2{iSong2}.covariance);
                end
            end

            cache1.mean1 = mean1;
            cache2.invCov2 = invCov2;

            if isfield(cache2, 'mean2')
                mean2 = cache2.mean2;
            else
                mean2 = zeros([nDim nSongs2]);
                for iSong2=1:nSongs2
                    mean2(:,iSong2) = features2{iSong2}.mean;
                end
                cache2.mean2 = mean2;
            end

            if twoSets
                if isfield(cache1, 'invCov1')
                    invCov1 = cache1.invCov1;
                else
                    invCov1 = zeros([nDim nDim nSongs1]);
                    for iSong1=1:nSongs1
                        invCov1(:,:,iSong1) = ...
                            inv(features1{iSong1}.covariance);
                    end
                end
                cache1.invCov1 = invCov1;

                % Compute non-symmetric KL-distance if two sets of feature
                % vectors were passed as arguments
                intermediateResult1 = zeros(nSongs1,1);
                intermediateResult2 = zeros([nSongs1 nDim]);
                for iSong=1:nSongs1
                    intermediateResult2(iSong,:) = mean1(:,iSong).'*invCov1(:,:,iSong);
                    intermediateResult1(iSong) = intermediateResult2(iSong,:) ...
                        * mean1(:,iSong);
                end

                distancematrix2 = zeros(nSongs1, nSongs2);
                for iSong=1:nSongs2
                    if ~silent, fprintf('Song %d of %d\n', iSong, nSongs2); end
                    m2sqr = mean2(:, iSong)*mean2(:, iSong)';
                    distancematrix2(:, iSong) = ...
                        reshape(sum(sum(invCov1.*features2{iSong}.covariance(:,:,ones(nSongs1,1)))), [], 1) + ...
                        reshape(sum(sum(invCov1.*m2sqr(:,:,ones(nSongs1,1)))), [], 1) + ...
                        intermediateResult1 - ...
                        intermediateResult2*(2*mean2(:, iSong)) - ...
                        nDim;
                end
                distancematrix2 = .5*distancematrix2;
            end

            % Compute non-symmetric KL-distance
            intermediateResult1 = zeros(nSongs2,1);
            intermediateResult2 = zeros([nSongs2 nDim]);

            for iSong=1:nSongs2
                intermediateResult2(iSong,:) = mean2(:,iSong).'*invCov2(:,:,iSong);
                intermediateResult1(iSong) = intermediateResult2(iSong,:) ...
                    * mean2(:,iSong);
            end

            distancematrix = zeros(nSongs2, nSongs1);
            
            for iSong=1:nSongs1
                if ~silent, fprintf('Song %d of %d\n', iSong, nSongs1); end
                m1sqr = mean1(:, iSong)*mean1(:, iSong)';
                distancematrix(:, iSong) = ...
                    reshape(sum(sum(invCov2.*features1{iSong}.covariance(:,:,ones(nSongs2,1)))), [], 1) + ...
                    reshape(sum(sum(invCov2.*m1sqr(:,:,ones(nSongs2,1)))), [], 1) + ...
                    intermediateResult1 - ...
                    intermediateResult2*(2*mean1(:, iSong)) - ...
                    nDim;
                
            end

            distancematrix = .5*distancematrix;
        end

        % Add non-symmetric distances to get the symmetric distance
        if twoSets
            distancematrix = distancematrix.' + distancematrix2;
        else
            distancematrix = distancematrix.' + distancematrix;
        end
        

        
      otherwise
        error(['Invalid type-argument '''  type '''.'])
    end


end
