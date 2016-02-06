%ISP_KMEANS  Binary tree k-means algorithm
%
% SYNTAX
%   [centroids, index]=isp_kmeans(vectors, nMeans)
%
% DESCRIPTION
%   ISP_KMEANS clusters a set of input vectors in 'nMeans' clusters.
%
% INPUT
%   vectors:
%     An MxN matrix of column vectors that are to be clustered.
%   nMeans:
%     The number of clusters the vectors are to be grouped into.
%
% OUTPUT
%   centroids:
%     An MxnMeans matrix where the column vectors are the centroids of
%     the clusters.
%   index:
%     The index of the cluster that each input vector belongs to.
%
% SEE ALSO
%   isp_gmmtrain.
%
% HISTORY
%   Created by Jesper Højvang Jensen in 2005.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [centroids, index]=isp_kmeans(vectors, nMeans)
   
    nVectors=size(vectors,2);
    centroids=zeros(size(vectors,1), nMeans);
    cluster=cell(1, nMeans);
    nVectorsInCluster=zeros(nMeans,1);
    centroids(:,1)=mean(vectors,2);
    nVectorsInCluster(1)=nVectors;
    cluster{1}=1:nVectorsInCluster(1);
    
    for n=2:nMeans
        [dummy, splitCluster]=max(nVectorsInCluster(1:n-1));
        vecIndices=cluster{splitCluster};
        nVecsInCluster=length(vecIndices);

        variance=sum(vectors(:, vecIndices).^2,2)/nVecsInCluster - centroids(:,splitCluster).^2;

        variance=sqrt(variance);
        centroids(:,n)=centroids(:,splitCluster)+variance;
        centroids(:,splitCluster)=centroids(:,splitCluster)-variance;


        for m=1:2
            t=(2*(centroids(:,splitCluster)-centroids(:,n)))'*vectors(:,vecIndices) > ...
                sum(centroids(:,splitCluster).^2) - sum(centroids(:,n).^2);
            % The statement above is equivalent to the commented one below, although
            % the former might be more tricky to understand
            % t=sum((vectors(:,vecIndices)-centroids(:,splitCluster*ones(1,nVecsInCluster))).^2,1) ...
            %     < sum((vectors(:,vecIndices)-centroids(:,n*ones(1,nVecsInCluster))).^2,1);

            cluster{n}=vecIndices(~t);
            cluster{splitCluster}=vecIndices(t);
            
            centroids(:,n)=mean(vectors(:, cluster{n}),2);
            centroids(:,splitCluster)=mean(vectors(:, cluster{splitCluster}),2);
        end
        
        nVectorsInCluster(n)=length(cluster{n});
        nVectorsInCluster(splitCluster)=length(cluster{splitCluster});
    end
    
    index=zeros(nVectors,1);
    for n=1:nMeans
        index(cluster{n})=n;
    end

end
