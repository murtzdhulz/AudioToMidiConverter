% ISP_EVALUATE  Evaluate a musical distance measures
%
% SYNTAX
%   [results, options] = isp_evaluate(distancemeasure, opt, 'field', 'value' ...)
%
% DESCRIPTION
%   Evaluates a distance measure on different test
%   collections. Possible collections include synthesized MIDI files,
%   the ISMIR 2004 genre classification training set, and the artist20 and
%   covers80 sets by Dan Ellis. It is necessary to download the test
%   collections and/or sondfonts for MIDI synthesis separately.
%
% INPUT
%   distancemeasure:
%     Either one struct or a cell array of structs specifying the
%     distance measures to be evaluated.
%   opt, field/value pairs:
%     A number of options either specified as fields of the struct 'opt'
%     or as "'field', value" pairs. Possible fields are:
%     experiment:
%       String or cell array of strings specifying which tests to
%       perform. Possible values are:
%       'artist20':
%         Artist recognition using Dan Ellis' artist20 data set.
%       'ballroom':
%         Style classification using the ISMIR 2004 ballroom set.
%       'bandwidth':
%         Measure sensitivity of instrument and melody recognition
%         accuracy for synthesized MIDI files to different, mixed
%         bandwidths. 
%       'bitrate':
%         Measure sensitivity of instrument and melody recognition
%         accuracy to MP3 compression for synthesized MIDI files. 
%       'covers80':
%         Cover identification using Dan Ellis' covers80 data set.
%       'downsample':
%         Measure sensitivity of instrument and melody recognition
%         accuracy to reduced bandwidth for synthesized MIDI files. 
%       'duration':
%         Measure sensitivity of instrument and melody recognition
%         accuracy to different tempi for synthesized MIDI files. 
%       'instrumentmelody':
%         Measure sensitivity of instrument and melody recognition
%         accuracy for synthesized MIDI files. 
%       'ismirtrainingset':
%         Genre classification using the ISMIR 2004 training set.
%       'multipleinstrument':
%         Measure sensitivity of instrument and melody recognition
%         accuracy for synthesized MIDI files when several instruments
%         play together. 
%       'silence':
%         Measure sensitivity of instrument and melody recognition
%         accuracy to added silence for synthesized MIDI files. 
%       'snr':
%         Measure sensitivity of instrument and melody recognition
%         accuracy to additive noise for synthesized MIDI files. 
%       'transpose':
%         Measure sensitivity of instrument and melody recognition
%         accuracy to transpositions for synthesized MIDI files. 
%       'all':
%         Perform all tests.
%       'allmidi'
%         Perform all MIDI based tests.
%
%     midiset:
%       Set to either 'long' or 'short', specifying whether to use MIDI
%       files of 3 minutes or 30 s length.
%     nInstruments:
%       Number of different instruments to use in the MIDI tests.
%       Default: 30;
%     nMidifiles = inf;
%       Number of MIDI files to use in the MIDI tests.
%       Default: 30;
%     percussion:
%       Boolean specifying whether to retain percussion in the MIDI
%       files. Default: false, i.e., remove percussion.
%     addsilence:
%       Vector specifying the amounts of silence added to songs for the
%       silence test. The amounts are specified as fractions of the
%       original song length. Default: [0 .05 .1 .15 .2].
%     duration:
%       Durations of songs for the duration test. Values are fractions of
%       the original song lengths. Default: [.5 .75 .9 1 1.1 1.25 1.5].
%     bandwidth:
%       Bandwidths for the bandwidth test.
%       Defaults: 0.5*[8000 11025 16000 22050 44100].
%     bitrate:
%       Bitrates for the bitrate test. Defaults: [32 40 48 56 64 inf].
%     snr:
%       SNRs for the snr test. Defaults: [-5 0 5 10 15 20 inf].
%     soundfont:
%       Cell array of strings specifying the sound fonts to use for the
%       MIDI tests. Default: 'isp_toolboxpath'/FluidR3 GM.SF2.
%     soundfontLabel:
%       Cell array of strings where the n'th element is a short text
%       describing the n'th sound font in 'soundfont'. The short texts
%       are used in legends. Default: 'Fluid'.
%     transpose:
%       Numbers of semitones to transpose for the transpose test.
%       Default: [-24 -19 -14 -10 -5 0 5 10 14 19 24].
%
%     artist20path, covers80path:
%       Paths to Dan Ellis' artist20 and covers80 data sets.
%       Defaults: './artist20' and './covers80', respectively.
%     ismir2004ballroompath, ismirgenrepath:
%       Paths to ISMIR 2004 ballroom data set and genre classification
%       training set. Defaults: '.' and './ismirgenre', respectively.
%     dataPath:
%       Path where results and optionally job files are stored.
%       Default: 'evaluationdata'.
%
%     distribute:
%       Either a number telling what should be distributed or 'ask'. If
%       set to 0, everything is executed immediately; if 1, distribute
%       feature extraction; if 2, merge features and distribute distance
%       computations; if 3, merge features and compute distances; if 4,
%       merge distances and extract results; if 5, show results.
%       Default: 'ask'.
%     distributeFunctions:
%       If MCC is used to compile the MATLAB code when distributing
%       stuff, and it cannot identify all functions, e.g. since they are
%       called using feval, additional function names can be specified as
%       a cell array of strings.
%     distributedDistancesPath:
%       Path where working files are stored when distributing distance
%       computations. Default: 'dataPath'/distributeddistances
%     distributedFeaturesPath:
%       Path where working files are stored when distributing feature
%       computations. Default: 'dataPath'/distributedfeatures
%     featuresPerJob:
%       The number of feature computations per job when distributing
%       feature computation. Default: 'auto'.
%     initializeDistributeDistances, initializeDistributeFeatures:
%       Booleans specifying whether isp_distribute shall be initialized
%       by isp_evaluate. If set to false, the user should have done this
%       manually before calling isp_evaluate. You probably wont need this
%       option. Defaults: false.
%
%     savedistancematrix:
%       Boolean specifying whether computed distance matrices should be saved.
%       Default: false.
%     savedistancematrixpath:
%       Path where computed distance matrices are stored. Default: 'dataPath'.
%     savefeature:
%       Boolean specifying whether extracted features should be saved.
%       Default: false.
%     savefeaturepath:
%       Path where extracted features are stored. Default: 'dataPath'.
%     savemidi:
%       Boolean specifying whether to save generated MIDI files, e.g. for
%       listening to them. Default: False.
%     savemidipath:
%       String specifying the path where MIDI files are saved. 
%       Default: 'savedmidi' under the data path.
%     savewav, savewavpath:
%       Similar to savemidi and savemidipath, just for the generated wav files.
%     usesavedwav:
%       If wav files have previously been saved, speed things up by
%       using the saved versions instead of re-generating them.
%       Default: false;
%
% OUTPUT
%   results:
%     A cell array where element (i, j) is a struct with evaluation
%     results of distance measure i evaluated in experiment j.
%   options:
%     A cell array of structs specifying the test options associated with
%     a result. It is basically the specified options supplemented with
%     default values.
%
% EXAMPLE
%   Evaluate the performance of the distance measure specified in the
%   isp_tichroma function:
%     isp_evaluate(isp_tichroma, 'experiment', 'covers80')
%
%   Also see isp_evaluatedemo1 and isp_evaluatedemo2 for examples of how
%   to use isp_evaluate.
%
% SEE ALSO
%   isp_evaluatedemo1, isp_evaluatedemo2.
%
% HISTORY
%   Created by Jesper H. Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.



function [resultArray, optionArray] = isp_evaluate(distancemeasure, varargin)

    %%
    %% Define constants and default values
    %%

    EXECUTE = 0;
    DISTRIBUTEFEATURES = 1;
    DISTRIBUTEDISTANCES = 2;
    COMPUTEDISTANCES = 3;
    EXTRACTRESULTS = 4;
    SHOWRESULTS = 5;

    INSTRUMENTMELODYTEST = 'instrumentmelody';
    TRANSPOSITIONTEST = 'transpose';
    BITRATETEST = 'bitrate';
    BANDWIDTHTEST = 'bandwidth';
    DOWNSAMPLETEST = 'downsample';
    SNRTEST = 'snr';
    SILENCETEST = 'silence';
    DURATIONTEST = 'duration';
    MULTIPLEINSTRUMENTTEST = 'multipleinstrument';
    ISMIRTRAININGSETTEST = 'ismirtrainingset';
    ARTIST20TEST = 'artist20';
    COVERS80TEST = 'covers80';
    BALLROOMTEST = 'ballroom';

    % Sound font settings
    testoptions.soundfont = {fullfile(isp_toolboxpath, 'FluidR3 GM.SF2')};
    testoptions.soundfontLabel = {'Fluid'};

    % Test settings
    testoptions.percussion = false;
    testoptions.savemidi = false;
    testoptions.savemidipath = '';
    testoptions.savewav = false;
    testoptions.savewavpath = '';
    testoptions.usesavedwav = false;
    testoptions.savefeature = false;
    testoptions.savefeaturepath = '';
    testoptions.savedistancematrix = false;
    testoptions.savedistancematrixpath = '';
    
    testoptions.transpose = round(linspace(-24,24,11));
    testoptions.bandwidth = 0.5*[8000 11025 16000 22050 44100];
    testoptions.bitrate = [32 40 48 56 64 inf];
    testoptions.snr = [-5 0 5 10 15 20 inf];
    testoptions.addsilence = [0 .05 .1 .15 .2];
    testoptions.duration = [.5 .75 .9 1 1.1 1.25 1.5];
    testoptions.nGenreFiles = inf;
    testoptions.nCoverFiles = inf;
    testoptions.maxSamples = inf;
    testoptions.removeSilence = true;
    testoptions.artist20path = './artist20';
    testoptions.covers80path = './covers80';
    testoptions.ismirgenrepath = './ismirgenre';
    testoptions.ismir2004ballroompath = '.';
    testoptions.nArtistFiles = inf;
    testoptions.nInstruments = inf;
    testoptions.nMidifiles = inf;
    testoptions.experiment = 'all';
    testoptions.midiset = 'short';

    % Settings for distributed computing
    testoptions.featuresPerJob = 'auto';
    testoptions.distributeFunctions = {};
    testoptions.distribute = 'ask'; % 'ask': Ask user
                                    % 0: Don't distribute anything
                                    % 1: Distribute feature extraction
                                    % 2: Merge features and distribute
                                    %    distance computation
                                    % 3: Merge features and compute distances
                                    % 4: Merge distances and extract results
                                    % 5: Show results

    testoptions.initializeDistributeDistances = true;
    testoptions.initializeDistributeFeatures = true;
    testoptions.distributedFeaturesPath = '';
    testoptions.distributedDistancesPath = '';
    testoptions.dataPath = fullfile('evaluationdata', '');

    %%
    %% Interpret input arguments
    %%

    if ~iscell(distancemeasure)
        distancemeasure = {distancemeasure};
    end
    nDistancemeasures = length(distancemeasure);

    testoptions = isp_interpretarguments(testoptions, varargin{:});

    % Set more defaults
    testoptions.distributedFeaturesPath = isp_nonempty(...
        testoptions.distributedFeaturesPath, ...
        fullfile(testoptions.dataPath, 'distributedfeatures', ''));
    testoptions.distributedDistancesPath = isp_nonempty(...
        testoptions.distributedDistancesPath, ...
        fullfile(testoptions.dataPath, 'distributeddistances', ''));
    testoptions.savefeaturepath = isp_nonempty(testoptions.savefeaturepath, ...
                                               testoptions.dataPath);
    testoptions.savedistancematrixpath = isp_nonempty(...
        testoptions.savedistancematrixpath, testoptions.dataPath);
    testoptions.savewavpath = isp_nonempty(testoptions.savewavpath, ...
                               fullfile(testoptions.dataPath, 'savedwav', ''));
    testoptions.savemidipath = isp_nonempty(testoptions.savemidipath, ...
            fullfile(testoptions.dataPath, 'savedmidi', ''));

    if numel(testoptions.soundfont) > numel(testoptions.soundfontLabel)
        warning(['More soundfonts specified than soundfont labels. Using ' ...
                 'soundfont file names as labels.'])
        for n=numel(testoptions.soundfontLabel)+1:numel(testoptions.soundfont)
            [dummy, testoptions.soundfontLabel{n}] = ...
                fileparts(testoptions.soundfont{n});
        end
    end

    if numel(testoptions.soundfontLabel) > numel(testoptions.soundfont) 
        testoptions.soundfontLabel(numel(testoptions.soundfont)+1:end) = [];
    end

    if ischar(testoptions.soundfont)
        testoptions.soundfont = {testoptions.soundfont};
    end

    allResultsFile = fullfile(testoptions.dataPath, 'all_results.mat');

    %
    % Initialize job distribution
    %
    if isequal(testoptions.distribute, 'ask') 
        % Ask what to do
        fprintf(1, 'Options:\n')
        fprintf(1, '0. Execute without distributing anything.\n')
        fprintf(1, '1. Distribute feature calculation.\n')
        fprintf(1, '2. Merge features and distribute distances.\n')
        fprintf(1, '3. Merge features, compute distances and extract results.\n')
        fprintf(1, '4. Merge distances and extract results.\n')
        fprintf(1, '5. Show results.\n')
        answer = str2num(input('What do you want to do? ', 's'));
        testoptions.distribute = answer;
    end

    switch testoptions.distribute
      case EXECUTE
      case DISTRIBUTEFEATURES
        isp_distribute('state', 1)
        if testoptions.initializeDistributeFeatures
            isp_distribute('start', testoptions.distributedFeaturesPath);
        end
        warning('off', 'ISP:noDistanceMatrix')
      case DISTRIBUTEDISTANCES
        isp_distribute('state', 1)
        if testoptions.initializeDistributeFeatures
            isp_distribute('read', testoptions.distributedFeaturesPath);
        end
        isp_distribute('state', 2)
        if testoptions.initializeDistributeDistances
            isp_distribute('start', testoptions.distributedDistancesPath);
        end
        warning('off', 'ISP:noDistanceMatrix')
      case COMPUTEDISTANCES
        isp_distribute('state', 1)
        if testoptions.initializeDistributeFeatures
            isp_distribute('read', testoptions.distributedFeaturesPath);
        end
      case EXTRACTRESULTS
        isp_distribute('state', 1)
        if testoptions.initializeDistributeFeatures
            isp_distribute('read', testoptions.distributedFeaturesPath);
        end
        isp_distribute('state', 2)
        if testoptions.initializeDistributeDistances
            isp_distribute('read', testoptions.distributedDistancesPath);
        end
      case SHOWRESULTS
        % Just show the results and quit
        load(allResultsFile)
        isp_plotresults(resultArray)
        return
      otherwise
        error('Invalid distribute command')
    end

    for iDistanceMeasure = 1:nDistancemeasures
        if isfield(distancemeasure{iDistanceMeasure}, 'usedFunctions')
            functionNames = distancemeasure{iDistanceMeasure}.usedFunctions;
            testoptions.distributeFunctions(end+1:end+length(functionNames)) ...
                = functionNames(:);
        end
    end
    testoptions.distributeFunctions = unique(testoptions.distributeFunctions);

    if ~exist(testoptions.dataPath, 'dir'), mkdir(testoptions.dataPath); end

    % Define test type
    if isequal(testoptions.experiment, 'all')
        experiment = {INSTRUMENTMELODYTEST, TRANSPOSITIONTEST, BITRATETEST, ...
                      BANDWIDTHTEST, DOWNSAMPLETEST, SNRTEST, SILENCETEST, ...
                      DURATIONTEST, MULTIPLEINSTRUMENTTEST, ISMIRTRAININGSETTEST, ...
                      ARTIST20TEST, COVERS80TEST, BALLROOMTEST};
    elseif isequal(testoptions.experiment, 'allmidi')
        experiment = {INSTRUMENTMELODYTEST, TRANSPOSITIONTEST, BITRATETEST, ...
                      BANDWIDTHTEST, DOWNSAMPLETEST, SNRTEST, SILENCETEST, ...
                      DURATIONTEST, MULTIPLEINSTRUMENTTEST};
    elseif isstr(testoptions.experiment)
        experiment = {testoptions.experiment};
    else
        experiment = testoptions.experiment;
    end

    % DOWNSAMPLETEST and BANDWIDTHTEST uses the same features. If both
    % are to be performed, the following block ensures they are executed
    % immediately after each other in order to reuse the features.
    if any(strcmp(experiment, BANDWIDTHTEST)) ...
            && any(strcmp(experiment, DOWNSAMPLETEST))
        performBothBandwidthTests = true;
        downsampleMask = strcmp(experiment, DOWNSAMPLETEST);
        experiment(downsampleMask) = [];
        bandwidthIdx = find(strcmp(experiment, BANDWIDTHTEST));
        experiment(bandwidthIdx+1:end+1) = experiment(bandwidthIdx:end);
        experiment{bandwidthIdx+1} = DOWNSAMPLETEST;
    else
        performBothBandwidthTests = false;
    end

    testoptions.experiment = experiment;

    if testoptions.savemidi, mkdir(testoptions.savemidipath); end
    if testoptions.savewav, mkdir(testoptions.savewavpath); end

    %%
    %% Perform actual test
    %%

    resultArray = {}; % For storing results
    optionArray = {};
    
    % Outer loop over the distance measures
    for iDistanceMeasure = 1:nDistancemeasures
        % Ensure all experiments are identical
        rand('state', 0);
        randn('state', 0);
        
        fprintf(1, '\nTesting distance measure %s.\n', ...
                distancemeasure{iDistanceMeasure}.name);
        clear feature feature0 features distancematrix0 distancematrices
            
        % Inner loop over the test types
        for iTest=1:numel(experiment(:))
            testtype = experiment{iTest};
            fprintf(1, 'Performing %s test.\n', testtype);

            options.distancemeasure = distancemeasure{iDistanceMeasure};


            resultFile = fullfile(testoptions.dataPath, ...
                [options.distancemeasure.name '_' testtype '_result.mat']);

            % Remove bandwidths larger than half the input sample rate
            testoptions.bandwidth( ...
                testoptions.bandwidth > 0.5*options.distancemeasure.samplerate ...
                ) = [];

            switch testtype
              case INSTRUMENTMELODYTEST
                [result, options] = isp_testskeleton(testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'midiInstruments', 'single', ...
                    'evalType', 'mixsets', ...
                    'modParameter', 'soundfont', ...
                    'modReference', '', ...
                    'testname', testtype);
                [result.modValues] = deal(testoptions.soundfontLabel);

              case MULTIPLEINSTRUMENTTEST
                [result, options] = isp_testskeleton(...
                    testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'iterator', 'soundfont', ...
                    'midiInstruments', 'multiple', ...
                    'evalType', 'overlap', ...
                    'testname', testtype);
                [result.iterationValue] = deal(testoptions.soundfontLabel{:});

              case TRANSPOSITIONTEST
                [result, options] = isp_testskeleton(...
                    testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'iterator', 'soundfont', ...
                    'midiInstruments', 'single', ...
                    'evalType', 'reference', ...
                    'modParameter', 'transpose', ...
                    'modReference', 0, ...
                    'testname', testtype);
                [result.iterationValue] = deal(testoptions.soundfontLabel{:});

              case BITRATETEST
                [result, options] = isp_testskeleton(...
                    testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'iterator', 'soundfont', ...
                    'midiInstruments', 'single', ...
                    'evalType', 'reference', ...
                    'modParameter', 'bitrate', ...
                    'modReference', inf, ...
                    'testname', testtype);
                [result.iterationValue] = deal(testoptions.soundfontLabel{:});

              case BANDWIDTHTEST
                [result, options] = isp_testskeleton(...
                    testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'iterator', 'soundfont', ...
                    'midiInstruments', 'single', ...
                    'evalType', 'reference', ...
                    'modParameter', 'bandwidth', ...
                    'modReference', 44100, ...
                    'testname', testtype);
                [result.iterationValue] = deal(testoptions.soundfontLabel{:});

              case DOWNSAMPLETEST
                [result, options] = isp_testskeleton(...
                    testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'iterator', 'soundfont', ...
                    'midiInstruments', 'single', ...
                    'evalType', 'nomix', ...
                    'modParameter', 'bandwidth', ...
                    'testname', testtype);
                [result.iterationValue] = deal(testoptions.soundfontLabel{:});

              case SNRTEST
                [result, options] = isp_testskeleton(...
                    testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'iterator', 'soundfont', ...
                    'midiInstruments', 'single', ...
                    'evalType', 'reference', ...
                    'modParameter', 'snr', ...
                    'modReference', inf, ...
                    'testname', testtype);
                [result.iterationValue] = deal(testoptions.soundfontLabel{:});

              case SILENCETEST
                [result, options] = isp_testskeleton(...
                    testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'iterator', 'soundfont', ...
                    'midiInstruments', 'single', ...
                    'evalType', 'reference', ...
                    'modParameter', 'addsilence', ...
                    'modReference', 0, ...
                    'testname', testtype);
                [result.iterationValue] = deal(testoptions.soundfontLabel{:});

              case DURATIONTEST
                [result, options] = isp_testskeleton(...
                    testoptions, ...
                    'dataset', 'midi', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'instrument', 'melody'}, ...
                    'iterator', 'soundfont', ...
                    'midiInstruments', 'single', ...
                    'evalType', 'reference', ...
                    'modParameter', 'duration', ...
                    'modReference', 1, ...
                    'testname', testtype);
                [result.iterationValue] = deal(testoptions.soundfontLabel{:});

              case COVERS80TEST
                [result, options] = isp_testskeleton(testoptions, ...
                    'dataset', 'covers80', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'title', 'title'}, ...
                    'evalType', 'nomix', ...
                    'testname', testtype, ...
                     'testset', {'all', 1}, ...
                     'trainingset', {'minusone', 2}, ...
                     'testsetname', {'title', 'split'});

              case BALLROOMTEST

                [result, options] = isp_testskeleton(testoptions, ...
                    'dataset', 'ismir2004ballroom', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'style', 'style', 'style', 'style'}, ...
                    'evalType', 'nomix', ...
                    'testname', testtype, ...
                     'testset', {'all', 'test', 'all', 'test'}, ...
                     'trainingset', {'minusone', 'training', 'minusone', 'training'}, ...
                     'testsetname', {'rhythm', 'split', 'difftempo', 'difftemposplit'}, ...
                     'hook', {[], [], @ballroomhook, @ballroomhook});

              case ISMIRTRAININGSETTEST
                [result, options] = isp_testskeleton(testoptions, ...
                    'dataset', 'ismirgenre', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalType', 'nomix', ...
                    'evalParameter', {'genre', 'genre'}, ...
                    'testset', {'all', 'all'}, ...
                    'trainingset', {'minusone', 'minusartist'}, ...
                    'testsetname', {'genre', 'artistfilter'}, ...
                    'testname', testtype);

              case ARTIST20TEST
                [result, options] = isp_testskeleton(testoptions, ...
                    'dataset', 'artist20', ...
                    'evalType', 'nomix', ...
                    'distancemeasure', distancemeasure{iDistanceMeasure}, ...
                    'evalParameter', {'artist', 'artist', 'artist', 'artist', 'artist', 'artist', 'artist'}, ...
                    'testset', {'all', 'all', 1, 2, 3, 4, 5}, ...
                    'trainingset', {'minusone', 'minusalbum', {2,3,4,5}, {1,3,4,5}, {1,2,4,5}, {1,2,3,5}, {1,2,3,4}}, ...
                    'testsetname', {'artist', 'albumfilter', 'set1', 'set2', 'set3', 'set4', 'set5'}, ...
                    'testname', testtype);
                                                     

              otherwise
                error('Invalid test type')
            end
            resultArray{iDistanceMeasure, iTest} = result;
            optionArray{iDistanceMeasure, iTest} = options;


            % Save computations
            if any(options.distribute==[EXECUTE, COMPUTEDISTANCES, EXTRACTRESULTS])
                save(resultFile, 'result', 'options', '-V6')

            end
            

        end

    end

    if any(options.distribute==[EXECUTE, COMPUTEDISTANCES, EXTRACTRESULTS])
        save(allResultsFile, 'resultArray', 'optionArray', '-V6')
    end

end


function [dstMtx, songs1, songs2] = ballroomhook(dstMtx, songs1, songs2)
    errorMargin=1.04;
    tempo1 = cat(1, songs1.tempo);
    tempo2 = cat(1, songs2.tempo);
    tmp1 = repmat(tempo1, 1, length(songs2));
    tmp2 = repmat(tempo2', length(songs1), 1);
    similarTempo = tmp1/errorMargin <= tmp2 & ...
        tmp2 < tmp1*errorMargin;
    [styles1{1:length(songs1)}] = songs1.style;
    [styles2{1:length(songs2)}] = songs2.style;
    sameStyle = strcmp(repmat(styles1(:), 1, length(songs2)), ...
                       repmat(styles2(:)', length(songs1),1));
    dstMtx(similarTempo & sameStyle) = nan;
    noMatchesIdx = all(similarTempo | ~sameStyle, 2);
    dstMtx(noMatchesIdx, :) = [];
    songs1(noMatchesIdx) = [];
end
