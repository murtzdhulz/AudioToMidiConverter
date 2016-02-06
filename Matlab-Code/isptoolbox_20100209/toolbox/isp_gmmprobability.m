%ISP_GMMPROBABILITY  Probability of vectors according to a specified Gaussian Mixture Model.
%
% SYNTAX
%   prob=isp_gmmprobability(gmmstruct, vectors)
%
% DESCRIPTION
%   Return the probability of vectors assuming the specified Gaussian
%   mixture model.
%
% INPUT
%   gmmstruct:
%     Struct with the following fields:
%     weight:
%       Vector of probabilities of each Gaussian mixture.
%     mean:
%       Matrix of column vectors containing the means of the individual
%       multivariate Gaussian components.
%     covariance:
%       3-dimensional matrix which contain the correlation matrices of
%       the individual multivariate Gaussian components.
%   vectors:
%     Cell array containing column vectors or matrix of column vectors.
%
% OUTPUT
%   prob:
%     Vector of probabilities of the vectors specified in 'vectors'
%
% SEE ALSO
%   isp_gmmrand, isp_gmmtrain.
%
% HISTORY
%   Created by Jesper H. Jensen 2005

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.

function prob=isp_gmmprobability(gmmstruct, vectors)
    
    % The number of mean value vectors to subtract at a time.
    blocksize=500;

    weight=gmmstruct.weight;
    meanVecs=gmmstruct.mean;
    diagonal=isfield(gmmstruct, 'diagonal') && gmmstruct.diagonal;
    
    nVectors=size(vectors,2);
    vecDim=size(vectors,1);
    nMixtures=gmmstruct.nMixtures;
    
    if size(weight,2) > size(weight,1)
        weight=weight';
    end

    if iscell(vectors)
        vectors=cell2mat(vectors);
    end

    
    if diagonal
        invDiagonals=1./gmmstruct.covariance;
        detCovariance=prod(gmmstruct.covariance, 1);
    else
        detCovariance=zeros(nMixtures,1);
        inverseCovariance=zeros(size(gmmstruct.covariance));
        for n=1:nMixtures
            detCovariance(n)=det(gmmstruct.covariance(:,:,n));
            inverseCovariance(:,:,n)=inv(gmmstruct.covariance(:,:,n));
        end
    end

    % Allocate variables
    prob=zeros(nVectors,1);
    exponents=zeros(1, nVectors);
    zeroMeanVectors=zeros(size(vectors));

    if diagonal
        vectorsSquared = vectors.^2;
    end

    % Calculate probabilities one mixture at a time
    for n=1:nMixtures
        % Ignore this mixture if the covariance matrix is singular
        if detCovariance(n)==0
            warning('Ignoring singular covariance matrix.')
            continue
        end
        constantFactor = 1./sqrt( ((2*pi)^vecDim)*detCovariance(n));

        if diagonal
            % Diagonal covariance
            d = -.5*invDiagonals(:,n);
            exponents = d'*vectorsSquared + d'*meanVecs(:,n).^2 - 2*(d.*meanVecs(:,n))'*vectors;
            prob = prob + (weight(n)*constantFactor)*exp(exponents');
        else
            % Full covariance
            if nVectors <= blocksize
                zeroMeanVectors = vectors - meanVecs(:,n*ones(nVectors,1));
            else
                tempVec=meanVecs(:,n*ones(blocksize,1));
                for m=1:blocksize:nVectors-blocksize;
                    zeroMeanVectors(:,m:m+blocksize-1)=vectors(:,m:m+blocksize-1)-tempVec;
                end
                zeroMeanVectors (:,m+blocksize:end)=vectors(:,m+blocksize:end)-meanVecs(:,n*ones(nVectors-m-blocksize+1,1));
            end
        
            exponents = -.5*sum(zeroMeanVectors.*(inverseCovariance(:,:,n)*zeroMeanVectors),1);
            prob = prob + (weight(n)*constantFactor)*exp(exponents');
        end
    end