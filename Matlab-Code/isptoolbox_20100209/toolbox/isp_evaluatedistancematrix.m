%isp_evaluatedistancematrix  Helper function for isp_evaluate
%
% SYNTAX
%   result = isp_evaluatedistancematrix(distanceMatrix, songs1, songs2, ...
%                                       classificationParameter, ignoreDiagonal)
%
% DESCRIPTION
%   This is a helper function for isp_evaluate. It extracts various
%   performance characteristics from a distance matrix and return them in
%   a struct.
%

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.



function result = isp_evaluatedistancematrix(distanceMatrix, ...
              songs1, songs2, classificationParameter, varargin)


    opts = isp_interpretarguments(struct('testset', '', 'trainingset', '', ...
                                         'hook', []), varargin{:});

    % distanceMatrix(i,j) is songs1(i) and songs2(j)
    nSongs1 = length(songs1);
    nSongs2 = length(songs2);

    if isempty(distanceMatrix)
        warning('ISP:noDistanceMatrix', 'Not evaluating distance matrix.')
        result=struct;
        return
    end

    fprintf(1, 'Evaluating distance matrix.\n')
    
    if ~isequal(size(distanceMatrix), [nSongs1 nSongs2])
        error('Mismatch between size of distance matrix and song arrays')
    end

    fractionOfFiniteDistances = mean(isfinite(distanceMatrix(:)));

    if ischar(opts.testset) && strcmp(opts.testset, 'all')
        d = distanceMatrix;
    else
        inSet=logical(zeros(nSongs1,1));
        if ~iscell(opts.testset)
            opts.testset={opts.testset};
        end
        for iSong=1:nSongs1
            for jSet=1:length(opts.testset);
                if isequal(songs1(iSong).dataset, opts.testset{jSet})
                    inSet(iSong) = true;
                end
            end
        end
        songs1(~inSet) = [];
        distanceMatrix(~inSet,:) = [];
        d = distanceMatrix;
    end


    if ischar(opts.testset)
        switch opts.trainingset
          case 'all'
            
          case 'minusone'
            fnames = intersect(fieldnames(songs1), ...
                               {'modification', 'soundfont'});
            if ~isequal(rmfield(songs1, fnames), rmfield(songs2, fnames))
                error('''minusone'' test set only valid when input song arrays are identical.')
            end
            d = distanceMatrix;
            d(1:nSongs1+1:end) = nan;
          otherwise
            if ischar(opts.trainingset) && ...
                    ~isempty(regexp(opts.trainingset, '^minus.+'))
                param=opts.trainingset(6:end);
                if ischar(songs1(1).(param))
                    % Speed this up
                    [paramVals1{1:nSongs1,1}] = songs1.(param);
                    [paramVals2{1:nSongs2,1}] = songs2.(param);
                    sameParam = strcmp(repmat(paramVals1, 1, nSongs2), repmat(paramVals2', nSongs1, 1));
                    distanceMatrix(sameParam) = nan;
                else
                    for iSong=1:nSongs2
                        for jSong=1:nSongs1
                            if isequal(songs2(iSong).(param), songs1(jSong).(param))
                                distanceMatrix(jSong, iSong) = nan;
                            end
                        end
                    end
                end
                d = distanceMatrix;
            end
        end
    else
        % Songs have a 'dataset' field
        inSet=logical(zeros(nSongs2,1));
        if ~iscell(opts.trainingset)
            opts.trainingset={opts.trainingset};
        end
        for iSong=1:nSongs2
            for jSet=1:length(opts.trainingset);
                if isequal(songs2(iSong).dataset, opts.trainingset{jSet})
                    inSet(iSong) = true;
                end
            end
        end
        songs2(~inSet) = [];
        distanceMatrix(:, ~inSet) = [];
        d = distanceMatrix;
    end

    if ~isempty(opts.hook)
        [d, songs1, songs2] = opts.hook(d, songs1, songs2);
    end

    nSongs1 = length(songs1);
    nSongs2 = length(songs2);

    % Compute the accuracy of a nearest neighbor classifier
    [minVal, minidx] = min(d, [], 2);

    % Use random ordering instead of the first if several minima exists
    minValuesMask = (d==repmat(minVal, 1, size(d, 2)));
    isnanMask = isnan(minVal);
    minValuesMask(isnanMask, :) = isnan(d(isnanMask,:));
    [dummy,minidx]=max(minValuesMask + .1*rand(size(minValuesMask)), [], 2);

    % Compare songs to their nearest neighbors
% $$$         [label{1:nSongs1}] = songs1(:).(classificationParameter);
% $$$         [nnLabel{1:nSongs1}] = songs2(minidx).(classificationParameter);
    nnAccuracy = 0;
    [label{1:nSongs1}] = songs1(:).(classificationParameter);
    [nnLabel{1:nSongs1}] = songs2(minidx).(classificationParameter);
    for iSong1 = 1:nSongs1
        nnAccuracy = nnAccuracy + isequal(label{iSong1}, nnLabel{iSong1});
    end
    nnAccuracy = nnAccuracy / nSongs1;

    %
    % Compute confusion matrix and histogram
    %
    if  ischar(label{1}) || isnumeric(label{1})
        [label1{1:nSongs1}] = songs1(:).(classificationParameter);
        [label2{1:nSongs2}] = songs2(:).(classificationParameter);
        if isnumeric(label1{1})
            nElements = numel(label1{1});
            tmp=zeros(nSongs1, nElements);
            for i=1:nSongs1
                tmp(i,:) = label1{i}(:)';
            end
            [label1, dummy, song1category] = unique(tmp, 'rows');
            tmp=zeros(nSongs2, nElements);
            for i=1:nSongs2
                tmp(i,:) = label2{i}(:)';
            end
            [label2, dummy, song2category] = unique(tmp, 'rows');
        else
            [label1, dummy, song1category] = unique(label1);
            [label2, dummy, song2category] = unique(label2);
        end
        
        confusionStatistics.label1 = label1;
        confusionStatistics.label2 = label2;

        confusionStatistics.matrix = zeros(length(label1), length(label2));

        for n=1:length(song1category)
            nnToCategory = song2category(minidx(song1category==n));
            for m=1:length(nnToCategory)
                confusionStatistics.matrix(n, nnToCategory(m)) = ...
                    confusionStatistics.matrix(n, nnToCategory(m)) + 1;
            end
        end

        [dummy, rankMtx] = sort(distanceMatrix, 2);
        for n=1:nSongs1
            rankMtx(n, rankMtx(n,:)) = 1:nSongs2;
        end
        rankHistogram = zeros(numel(songs2), 1);
        for s1=1:length(label1)
            for s2=1:length(label2)
                if isequal(label1(s1), label2(s2));
                    t=rankMtx(song1category==s1, song2category==s2);
                    for n=1:numel(t)
                        rankHistogram(t(n)) = rankHistogram(t(n)) + 1;
                    end
                end
            end
        end
    else
        confusionStatistics = [];
        rankHistogram = [];
    end    
    

    % overlapping number of instruments
    overlapStatistics = [];
    if isfield(songs1(1), 'instrument')
        nInstrumentsPerSong = numel(songs1(1).instrument);
    else
        nInstrumentsPerSong = nan;
    end

    if strcmp(classificationParameter, 'instrument') && nInstrumentsPerSong > 1
        % Compute overlap statistics
        if ~isequal(songs1, songs2)
            warning(['Implicit assumptions about test songs set 1 and 2 violated.'])
        end

        [song1midifiles{1:nSongs1}] = songs1.melody;
        [song2midifiles{1:nSongs2}] = songs2.melody;
        [uniqueSongs1, idx1, idx2] = unique(song1midifiles);
        
        for iSong=1:length(uniqueSongs1)
            sameMidifile(iSong, 1:nSongs2) = strcmp(uniqueSongs1{iSong}, song2midifiles);
        end
        sameMidifile = sameMidifile(idx2,:);

        d = distanceMatrix;
        d(sameMidifile) = nan;
        [minVal, minidx] = min(d, [], 2);
        
        % Use random ordering instead of the first if several minima exists
        minValuesMask = (d==repmat(minVal, 1, size(d, 2)));
        isnanMask = isnan(minVal);
        minValuesMask(isnanMask, :) = isnan(d(isnanMask,:));
        [dummy,minidx]=max(minValuesMask + .1*rand(size(minValuesMask)), [], 2);

        % Compare songs to their nearest neighbors
        nIdenticalInstruments = zeros(nSongs1,1);
        for iSong1 = 1:nSongs1
            label = songs1(iSong1).instrument;
            nnLabel = songs2(minidx(iSong1)).instrument;
            nIdenticalInstruments(iSong1) = sum(label==nnLabel);
        end

        for i=0:nInstrumentsPerSong
            nearestNeighborAccuracy(i+1) = mean(nIdenticalInstruments==i);
        end

        overlapStatistics.nearestNeighborAccuracy = nearestNeighborAccuracy;
        overlapStatistics.description = ['''nearestNeighborAccuracy(i+1)'' is the fraction of nearest neighbors with i identical instruments.'];

    end
    
    result.nnAccuracy = nnAccuracy;
    result.fractionOfFiniteDistances = fractionOfFiniteDistances;
    result.overlapStatistics = overlapStatistics;
    result.confusionStatistics = confusionStatistics;
    result.rankHistogram = rankHistogram;
end
