%ISP_GMMTRAIN  Train a Gaussian mixture model from data
%
% SYNTAX
%   [gmmstruct, options]=isp_gmmtrain(vectors, options ...)
%
% DESCRIPTION
%   Trains a Gaussian mixture model (GMM) from a set of training vectors.
%   
% INPUT
%   vectors:
%     MxN matrix consisting of N M-dimensional column vectors.
%   options ...:
%     The following parameters can be set as field names in structs, or
%     be specified as field/value pairs:
%     nMixtures:
%       Number of mixtures (default: 10)
%     diagonal:
%       Boolean specifying if covariance matrices are diagonal (default:
%       true)
%     nIterations:
%       number of iterations of EM algorithm. Must be an integer >= 0 or
%       the string 'auto' (default: 'auto')
%     daem:
%       Boolean specifying if the deterministic annealing EM-algorithm is
%       used. (default: false)
%     maxIterations:
%       Maximum number of iterations if 'nIterations' is set to 'auto'
%       (default: 100 for full covariance matrices and 1000 for diagonal)
%     convergenceThreshold:
%       Theshold used as stopping criterion when nIterations is set to
%       'auto' (default: 0.2)
%     daemStart:
%       DAEM start value in the range 0 to 1 (default: 0.5)
%     daemStop:
%       DAEM stop value  in the range 0 to 1 and larger than daemStart
%       (default: 1)
%     daemGain:
%       Multiplicative increase of daem factor (default: 1.2)
%     
% OUTPUT
%   gmmstruct:
%     Structure describing the trained GMM.
%
% SEE ALSO
%   isp_gmmprobability, isp_gmmrand.
%
% HISTORY
%   2005: Created by Jesper H. Jensen
%   2006-03-30: First version for the ISP toolbox.
%   2007-04-02: Corrected a bug that caused invalid output in very rare cases.
%   2007-11-12: Restructuring according to new ISP toolbox layout.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [gmmstruct, options]=isp_gmmtrain2(vectors, varargin)

    if nargin==0, error('No input parameters specified'); end
    
    %
    % Default values
    %

    options = struct('nMixtures', 10, ...
                     'diagonal', true, ...
                     'nIterations', 'auto', ...
                     'daem', false, ...
                     'maxIterations', [], ...
                     'daemStart', 0.5, ...
                     'daemStop', 1, ...
                     'daemGain', 1.2, ...
                     'silent', false, ...
                     'convergenceThreshold', 0.2);
    
    %
    % Interpret input parameters
    %

    options = isp_interpretarguments(options, varargin{:});

    nVectors=size(vectors,2);
    vecDim=size(vectors,1);
    nMixtures=options.nMixtures;
    diagonal=options.diagonal;
    silent=options.silent;

    if isempty(options.maxIterations)
        if diagonal
            options.maxIterations = 1000;
        else
            options.maxIterations = 100;
        end
    end

    if strcmp(options.nIterations, 'auto')
        nIterations=options.maxIterations;
        convergenceThreshold=options.convergenceThreshold;
    else
        nIterations=options.nIterations;
        convergenceThreshold=-inf;
    end

    if options.daemStart <=0 || options.daemStart >1 || options.daemGain <=1
        error('Illegal value of DAEM parameters')
    end

    if nIterations==0 || ~options.daem
        DAEMbetaValues=1;
    else
        n=1;
        beta=options.daemStart;
        while beta < options.daemStop
            DAEMbetaValues(n)=beta;
            beta=beta*options.daemGain;
            n=n+1;
        end
        DAEMbetaValues(n)=options.daemStop;
    end

    if nMixtures==1
        % We have an exact solution. No need to do anything iteratively
        nIterations=0;
    end

    if iscell(vectors)
        vectors=cell2mat(vectors);
    end

   
    %
    % Use kmeans to initialize GMM
    %
    
    if nMixtures==1
        means=mean(vectors,2);
        indices(1:size(vectors,2))=1;
    else
        [means, indices]=isp_kmeans(vectors, nMixtures);
    end

    weights=zeros(1, nMixtures);
    for n=1:nMixtures
        clusterVectors=vectors(:, indices==n);
        numPoints=size(clusterVectors,2);
        if diagonal
            covariances(:,n) = mean((clusterVectors).*conj(clusterVectors),2) ...
                - means(:,n).*conj(means(:,n));
        else
            covariances(:,:,n) = clusterVectors*clusterVectors'/numPoints - ...
                means(:,n)*means(:,n)';
        end
        weights(n)=numPoints/nVectors;
    end

    % Set "very small value" which is added to a nearly-singular
    % covariance matrix to make it non-singular
    if diagonal
        delta=0.01*mean(covariances, 2);
    else
        delta=diag(diag(0.01*mean(covariances, 3)));
    end

    
    %
    % Do the EM-algorithm
    %

    % Allocate variables
    probVecFromMixture=zeros(nVectors,nMixtures);
    probability=zeros(1, nIterations);
    meansSubtracted=zeros(vecDim, nVectors);
    diagonalValues=zeros(vecDim, 1);

    % Reset things
    lastwarn('')
    nSingularities=0;
    loglikelihood=[];
    iterationNumber=0;

    if diagonal
        vectorsSqr = vectors.^2;
    end


    % DAEM loop. If DAEM isn't used, it is effectively skipped.
    for beta=DAEMbetaValues;

        % EM iterations
        while true

            % Halt if this is taking too long
            if  iterationNumber >= nIterations
                if  nIterations > 0 && iterationNumber == nIterations
                    if ~silent
                        fprintf('Interrupting EM-algorithm at step %d\n', iterationNumber)
                    end
                end
                break;
            end

            iterationNumber=iterationNumber+1;

            if ~silent
                if convergenceThreshold==-inf
                    fprintf(1, 'Iteration %d (max. %d).\n', iterationNumber, nIterations);
                else
                    fprintf(1, 'Iteration %d of %d.\n', iterationNumber, nIterations);
                end
            end

            if diagonal
                % Diagonal covariance matrixes
                diagproduct = prod(covariances, 1);
                maxv=max(covariances, [], 1);
                condMask = (maxv==0) | (diagproduct==0);
                if any(condMask)
                    covariances(:,condMask) = covariances(:,condMask) + delta(:,ones(sum(condMask),1));
                    diagproduct = prod(covariances);
                else
                    rcond = min(covariances, [], 1)./maxv;
                    condMask = rcond < eps;
                    if any(condMask)
                        nSingularities = nSingularities+1;
                        covariances(:,condMask) = covariances(:,condMask) + delta(:,ones(sum(condMask),1));
                        diagproduct = prod(covariances);
                    end
                end

                invDiag = -0.5./covariances;
                constfactor = weights./((2*pi)^(vecDim/2).*sqrt(diagproduct));

                probVecFromMixture = constfactor(ones(nVectors, 1),:).*exp( ...
                    vectorsSqr'*invDiag + ...
                    repmat(sum((means.^2).*invDiag,1), nVectors, 1) - ...
                    vectors'*(2*means.*invDiag));

            else
                % Full covariance matrixes
                for k=1:nMixtures
                    meansSubtracted=vectors - means(:, k*ones(nVectors,1));
                    detCov = det(covariances(:,:,k));
                    constfactor=weights(k)/((2*pi)^(vecDim/2)*sqrt(detCov));
                    probVecFromMixture(:,k)=constfactor*exp(-0.5* sum(meansSubtracted.*(covariances(:,:,k)\meansSubtracted), 1) );
                    if ~isempty(lastwarn) || real(detCov) < 0 || ~all(isfinite(probVecFromMixture(:,k)))
                        covariances(:,:,k)=covariances(:,:,k) + delta;
                        constfactor=weights(k)/((2*pi)^(vecDim/2)*sqrt(det(covariances(:,:,k))));
                        probVecFromMixture(:,k)=constfactor*exp(-0.5* sum(meansSubtracted.*(covariances(:,:,k)\meansSubtracted), 1) );
                        lastwarn('')
                        nSingularities=nSingularities+1;
                    end
                end
            end

            probyi=sum(probVecFromMixture,2);

            if any(probyi==0)
                %warning('Numerical problems ...')
                tmp = (probyi==0);
                probVecFromMixture(tmp, :) = realmin;
                probyi(tmp) = nMixtures*realmin;
            end
            
            if options.daem
                % Deterministic Annealing EM (DAEM) algorithm
                DAEMprobVecFromMixture=probVecFromMixture.^beta;
                DAEMprobyi=sum(DAEMprobVecFromMixture,2);
                lambda_ik=DAEMprobVecFromMixture./repmat(DAEMprobyi,1,nMixtures);
            else
                % Ordinary EM algorithm
                lambda_ik=probVecFromMixture./repmat(probyi,1,nMixtures);
            end

            lambda_k=sum(lambda_ik,1);
            weights=lambda_k/nVectors;

            if diagonal
                % Diagonal covariance matrixes
                means = (vectors*lambda_ik)./lambda_k(ones(vecDim,1), :);
                temp = vectorsSqr*lambda_ik ...
                       + means.^2.*lambda_k(ones(vecDim,1),:) ...
                       - 2*(vectors*lambda_ik).* means;
                covariances = temp./lambda_k(ones(vecDim,1), :);
            else
                % Full covariance matrixes
                for k=1:nMixtures
                    means(:,k)=vectors*lambda_ik(:,k)/lambda_k(k);
                    meansSubtracted=vectors - means(:, k*ones(nVectors,1));

                    p(1, :)=lambda_ik(:,k);
                    temp = ( p(ones(vecDim,1), :) .* meansSubtracted)*meansSubtracted';
                    covariances(:,:,k)=temp/lambda_k(k);
                end
            end

            loglikelihood(iterationNumber)=sum(log(probyi));

            if iterationNumber>10 && var(loglikelihood(iterationNumber-9:iterationNumber)) < convergenceThreshold
                if beta==1
                    fprintf('Stopping EM-algorithm at step %d\n', iterationNumber)
                else
                    fprintf('Finished using beta value of %f\n', beta)
                end
                break;
            end

        end

    end

    gmmstruct=struct('type', 'gmm', ...
        'weight',weights, ...
        'mean', means, ...
        'covariance', covariances, ...
        'nMixtures', nMixtures, ...
        'diagonal', diagonal);

    if nSingularities > 0
        fprintf('Number of singularities: %d\n', nSingularities);
    end

end
