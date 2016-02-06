%ISP_GMMRAND  Generate random vectors from a Gaussian Mixture Model.
%
% SYNTAX
%   vectors=isp_gmmrand(gmmstruct, nVectors, enforceMix)
%
% DESCRIPTION
%   Generate random vectors from a Gaussian mixture model.
%
% INPUT
%   gmmstruct:
%     Struct with the following fields (as returned by isp_gmmtrain):
%     weight:
%       Vector of probabilities of each Gaussian mixture.
%     mean:
%       Matrix of column vectors containing the means of the individual
%       multivariate Gaussian components.
%     diagonal:
%       True if covariance matrices are diagonal. False otherwise.
%     covariance:
%       If 'diagonal' is true, a 2-dimensional matrix where covariance(:,k)
%       is the diagonal of the k'th covariance matrix. Otherwise it is a
%       3-dimensional matrix where covariance(:,:,k) is the covariance
%       matrix of the k'th multivariate Gaussian.
%   nVectors:
%     The number of random vectors to return (default: 1).
%   enforceMix:
%     If false (default), let the choice of Gaussian be random,
%     distributed according to 'weight'. If true, let the number of
%     vectors from each mixture be proportional to 'weight'.
%
% OUTPUT
%   vectors:
%     Matrix where the columns consist of random vectors.
%
% SEE ALSO
%   isp_gmmprobability, isp_gmmtrain.
%
% HISTORY
%   Created by Jesper H. Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function vectors=isp_gmmrand(gmmstruct, nVectors, enforceMix)
    weights=gmmstruct.weight;
    means=gmmstruct.mean;
    covariances=gmmstruct.covariance;
    diagonal=gmmstruct.diagonal;

    %
    % Interpret arguments
    %
    
    if ~exist('nVectors', 'var')
        nVectors=1;
    end

    if ~exist('enforceMix', 'var')
        enforceMix=false;
    end

    if size(weights,2) > size(weights,1)
        weights=weights';
    end

    
    %
    % Allocate memory
    %
    
    vecDim=size(means,1);
    vectors=zeros(vecDim,nVectors);

    
    if enforceMix
        thresholds=(0:nVectors)/nVectors;
    else
        thresholds=rand(1,nVectors);
    end
    cumulativeWeights=[0; cumsum(weights)];
    
    
    if diagonal
        covarianceSqrt=sqrt(covariances);
    else
        covarianceSqrt=repmat(nan,size(covariances));
    end
    
    for n=1:nVectors

        % Use a uniformly distributed random number to determine which of the
        % Gaussian mixtures to use
        thisThreshold=thresholds(n);

        % Do a binary search
        maxIndex=length(cumulativeWeights);
        minIndex=1;
        while maxIndex-minIndex ~= 1
            middleIndex=ceil(0.5*(maxIndex + minIndex));
            if cumulativeWeights(middleIndex) < thisThreshold
                minIndex=middleIndex;
            else
                maxIndex=middleIndex;
            end
        end
        mixtureNumber=maxIndex-1;

        if diagonal
            vectors(:,n)=covarianceSqrt(:,mixtureNumber).*randn(vecDim,1) ...
                + means(:,mixtureNumber);
        else
            % Only calculate the "square root" of a correlation matrix once.
            if isnan(covarianceSqrt(1,1,mixtureNumber))
                covarianceSqrt(:,:,mixtureNumber) = ...
                    sqrtm(covariances(:,:,mixtureNumber));
            end

            % Calculate a vector from a multivariate Gaussian distribution by
            % multiplying the "square root" of the correlation matrix by a vector of
            % independent Gaussians and add the mean.
            vectors(:,n)=covarianceSqrt(:,:,mixtureNumber)*randn(vecDim,1) ...
                + means(:,mixtureNumber);
        end
    end
