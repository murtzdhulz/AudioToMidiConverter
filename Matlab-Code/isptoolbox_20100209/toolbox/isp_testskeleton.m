%ISP_TESTSKELETON  Internal function providing a framework for isp_evaluate.
%
% SYNTAX:
%   [result, options] = isp_testskeleton(options ...)
%
% DESCRIPTION
%   This is a function providing a framework for isp_evaluate. It is not
%   supposed to be called directly.
%
% SEE ALSO
%   isp_evaluate.
%
% HISTORY
%   Created by Jesper Højvang Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [result, options] = isp_testskeleton(varargin)

    options = isp_interpretarguments(true, struct( ...
        'testset', '', ...
        'trainingset', '', ...
        'testsetname', '', ...
        'modParameter', '', ...
        'modReference', '', ...
        'hook', '', ...
        'evalParameter', '', ...
        'evalType', '', ...
        'iterator', ''), ...
        varargin{:});

    if ~iscell(options.evalParameter)
        options.evalParameter = {options.evalParameter};
    end

    options.testset = isp_nonempty(options.testset, ...
                                 repmat({'all'}, size(options.evalParameter)));
    options.hook = isp_nonempty(options.hook, ...
                                 repmat({[]}, size(options.evalParameter)));
    options.trainingset = isp_nonempty(options.trainingset, ...
                            repmat({'minusone'}, size(options.evalParameter)));

    options.testsetname = isp_nonempty(options.testsetname, ...
                                       options.evalParameter);

    hasReference = ~isempty(options.modReference);
    if hasReference && ~strcmp(options.evalType, 'reference')
        error(['Reference parameter specified, '...
               'but test type isn''t ''reference''.'])
    end

    isModified = ~isempty(options.modParameter);
    if isModified
        modValues = options.(options.modParameter);
        if ~iscell(modValues)
            modValues = num2cell(modValues);
        end
    else
        modValues = [];
    end

    doIterate = ~isempty(options.iterator);
    if doIterate
        iterationValues = options.(options.iterator); 
        if ~iscell(iterationValues)
            iterationValues = num2cell(iterationValues);
        end
    else
        iterationValues = {[]};
    end

    [unmodSongArray, options]=isp_makesonglist(options);
    dummy={struct};
    if ~isfield(unmodSongArray, 'modification')
        [unmodSongArray.modification] = dummy{ones(size(unmodSongArray))};
    end

    dist = options.distancemeasure;

    distTime = 0;
    nDists = 0;


    for iIter=1:length(iterationValues);
        songArray = modifySongArray(unmodSongArray, options.iterator, ...
                                    iterationValues(iIter));

        % Specify songs
        songArray=modifySongArray(songArray, options.modParameter, ...
                                  modValues);

        % Compute features
        features = cell([numel(songArray) numel(options.modParameter)]);
        [features, songArray]=skeleton_extractfeature(options, songArray, dist);

        % Compute reference features if needed
        if hasReference
            refNum = 0;
            for n=1:length(modValues)
                if isequal(options.modReference, modValues{n})
                    refNum=n;
                    break
                end
            end

            if refNum == 0
                refSongArray = modifySongArray(unmodSongArray, ...
                                     options.iterator, iterationValues(iIter));
                refSongArray = modifySongArray(refSongArray, ...
                                 options.modParameter, {options.modReference});
                [refFeatures, refSongArray] = skeleton_extractfeature(...
                    options, refSongArray, dist);
            else
                refSongArray=songArray(:, refNum);
                refFeatures = features(:, refNum);
            end
        end

        % Compute distances between features and evaluate them
        switch options.evalType
          case {'mixsets', 'nomix', 'overlap'}
            
            distancematrix=cell(size(songArray,2));
            for iParam=1:size(songArray,2)
                switch options.evalType
                  case 'mixsets'
                    jParamValues = 1:size(songArray,2);
                  case {'nomix', 'overlap'}
                    jParamValues = iParam;
                  otherwise
                    error('fsgfjsd')
                end
                for jParam=jParamValues
                    [distancematrix{iParam, jParam}, t] = ...
                        skeleton_computedistance(options, dist, ...
                                                 features(:,jParam), features(:,iParam));
                    for iEval=1:numel(options.evalParameter)
                        result(iIter).(options.testsetname{iEval})(iParam, jParam) ...
                            = isp_evaluatedistancematrix(...
                                distancematrix{iParam, jParam}, ...
                                songArray(:,iParam), songArray(:,jParam), ...
                                options.evalParameter{iEval}, ...
                             'testset', options.testset{iEval}, ...
                             'trainingset', options.trainingset{iEval}, ...
                             'hook', options.hook{iEval});
                    end
                end
            end

          case 'reference'
            distancematrix=cell(1, size(songArray,2));
            for n=1:size(songArray,2)
                [distancematrix{n}, t]=skeleton_computedistance(...
                    options, dist, features(:,n), refFeatures);
                for iEval=1:numel(options.evalParameter)
                    result(iIter).(options.testsetname{iEval})(n) ...
                          = isp_evaluatedistancematrix(distancematrix{n}, ...
                             songArray(:,n), refSongArray, ...
                             options.evalParameter{iEval},...
                             'testset', options.testset{iEval}, ...
                             'trainingset', options.trainingset{iEval}, ...
                             'hook', options.hook{iEval});
                end
            end

          otherwise
            error('Invalid or unspecified evalType.');
        end
    end

    [result.distancename] = deal(dist.name);
    [result.test] = deal(options.testname);
    [result.iterator] = deal(options.iterator);
    if doIterate
        [result.iterationValue] = deal(iterationValues{:});
    else
        [result.iterationValue] = deal([]);
    end
    
    [result.evalType] = deal(options.evalType);
    [result.evalParameter] = deal(options.evalParameter);
    [result.modParameter] = deal(options.modParameter);
    [result.modValues] = deal(modValues);
    [result.testsetname] = deal(options.testsetname);
    
    if options.savefeature
        featurefilename = fullfile(options.savefeaturepath, ...
                              [dist.name '_' options.testname '_feature.mat']);
        includedVariables = {};
        for variableName={'result', 'options', 'songArray', 'features', ...
                          'refFeatures', 'songArray'}
            if exist(variableName{1}, 'var') && ~isempty(variableName{1})
                includedVariables(end+1) = variableName;
            end
        end
        save(featurefilename, includedVariables{:})
    end
    
    if options.savedistancematrix
        distancefilename = fullfile(options.savedistancematrixpath, ...
                             [dist.name '_' options.testname '_distance.mat']);
        includedVariables = {};
        for variableName={'result', 'options', 'songArray', ...
                          'distancematrix', 'songArray'}
            if exist(variableName{1}, 'var') && ~isempty(variableName{1})
                includedVariables(end+1) = variableName;
            end
        end
        save(distancefilename, includedVariables{:})
    end
    
end


function [features, song]=skeleton_extractfeature(testoptions, song, distmsr)

    song = song';
    features=cell(size(song));
    
    switch testoptions.distribute
      case 0
        [features, song]=isp_extractfeature(song, distmsr);
      case {1,2,3}
        if isequal(testoptions.featuresPerJob, 'auto')
            testoptions.featuresPerjob = ceil(20/size(song, 1))*size(song, 1);
        end

        cmd='[features(i:endidx), song(i:endidx)]=isp_extractfeature(song(i:endidx), distmsr);';
        isp_distribute('state', 1)
        for i=1:testoptions.featuresPerjob:numel(song)
            endidx = min(i+testoptions.featuresPerjob-1, numel(song));
            isp_distribute(cmd, 2);
        end

      case 4
        %features = {};
      otherwise
        error('Invalid value of testoptions.distribute')
    end

    features = features';
    song = song';
end


function varargout=skeleton_computedistance(testoptions, dist, features1, ...
                                            features2)
    varargout={[], 0};
    cmd=['[varargout{1:nargout}]=isp_computedistance(dist, features1, ' ...
         'features2);'];
    switch testoptions.distribute
      case {0,3}
        fprintf(1, 'Computing distance matrix.\n')
        eval(cmd);
      case 1
      case {2,4}
        isp_distribute('state', 2)
        isp_distribute(cmd, nargout);
      otherwise
        error('Invalid value of testoptions.distribute')
    end
    
end


function songArray=modifySongArray(songArray, modParameter, modValues)
    if isempty(modParameter)
        return
    end

    songArray = repmat(songArray(:), 1, numel(modValues));
    for n=1:numel(modValues)
        for m=1:size(songArray,1)
            if strcmp(modParameter, 'soundfont')
                songArray(m, n).soundfont = modValues{n};
            else
                songArray(m, n).modification.(modParameter) ...
                    = modValues{n};
            end
        end
    end
end