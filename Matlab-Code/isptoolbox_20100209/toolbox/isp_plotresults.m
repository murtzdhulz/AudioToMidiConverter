%ISP_PLOTRESULTS  Plots results returned by isp_evaluate
%
% SYNTAX
%   [h, plotargs, options] = isp_plotresults(results, options ...)
%
% DESCRIPTION
%   Plots results returned by isp_evaluate.
%
% INPUT
%   results:
%     A result struct or a cell array of structs as returned by
%     isp_evaluate. If it is a cell array, element (m,n) must be the
%     results from the m'th distance measure and the n'th test (as
%     returned by isp_evaluate).
%   options:
%     Structs or field/value pairs specifying some or all of the
%     following properties:
%     evaluations:
%       Mask specifying which evaluations to plots. Typically,
%       [true false] will show instrument only, and [false true] will
%       show melody only. Default: Show all results.
%     linestyleoptions:
%       Line style options passed on to isp_linestyle.
%     integration:
%       How to integrate results when the evaluation of a MIDI test has
%       been repeated using different sound fonts. Possible values are
%       'mean', 'stddev', 'maxmin' and 'all', plotting the mean accuracy,
%       the mean plus/minus the standard deviation, the lowest/highest
%       accuracy and all results, respectively.
%     xaxis:
%       What to plot on the x-axis. Possible values are 'modification'
%       and 'distancemeasure'. Default: 'distancemeasure'.
%     modLabel:
%       Label for the x-axis when 'xaxis' is set to 'modification'.
%     distLabel:
%       Label for the x-axis when 'xaxis' is set to 'distancemeasure'.
%     plotTitle:
%       Title of plot.
%
% OUTPUT
%   h:
%     Handles for all plots and subplots.
%   plotargs:
%     Arguments passed to call isp_plot to make the plots.
%   options:
%     The input options supplemented with default values.
%
% SEE ALSO
%   isp_evaluate, isp_linestyle and isp_plot.
%
% HISTORY
%   Created by Jesper Højvang Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.



function [h, plotargs, options] = isp_plotresults(results, varargin)

    if ~iscell(results)
        results={results};
        comparison = false;
    end

    for n=1:length(results)
        test{n} = results{n}(1).test;
    end
    
    
    % Ought to verify the structure of 'results'
    for iTest = 1:size(results, 2);
        testDescription = results{1, iTest}(1).test;
        default = struct;
        switch testDescription

          case 'instrumentmelody'
          case 'transpose'
            %        tests = {'instrument', 'melody'};
% $$$             testNames = {'Effect of transpositions on instrument accuracy', ...
% $$$                          'Effect of transpositions on melody accuracy'};
% $$$             name=[upper(testDescription(1)) testDescription(2:end)];
% $$$             modLabel = 'Number of semitones transposed';
% $$$             modLegendLabel = 'Transposing %s semitones';
            default.modLabel = 'Transposition in Semitones';
            
          case 'bitrate'
            default.modLabel = 'Bitrate in kbps';
          case 'bandwidth'
            default.modLabel = 'Bandwidth in Hz';
          case 'downsample'
            default.modLabel = 'Bandwidth in Hz';
          case 'snr'
            default.modLabel = 'SNR in dB';
          case 'silence'
          case 'duration'
          case 'multipleinstrument'
          case 'ismirtrainingset'
          case 'artist20'
	    for jTest=1:size(results(:, iTest), 1)
                split.nnAccuracy = mean([results{jTest, iTest}.set1.nnAccuracy
				       results{jTest, iTest}.set2.nnAccuracy
				       results{jTest, iTest}.set3.nnAccuracy
				       results{jTest, iTest}.set4.nnAccuracy
				       results{jTest, iTest}.set5.nnAccuracy]);
		split.overlapStatistics = [];
		results{jTest, iTest}.split = split;
		results{jTest, iTest}.testsetname(3) = {'split'};
		results{jTest, iTest}.testsetname(4:end) = [];
% 		results{jTest, iTest} = rmfield(results{jTest, iTest}, ...
% 				{'set1', 'set2', 'set3', 'set4', 'set5'});
	    end
          case 'covers80'
          case 'ballroom'
          otherwise
            warning(['Unknown test type ' testDescription '.'])
        end
        
        [h{iTest}, plotargs{iTest}, options(iTest)] = plotresultscomp_helper(results(:, iTest), default, ...
                                              varargin{:});


    end
end % plotresults





function [h, plotargs, options] = plotresultscomp_helper(results, varargin)

    testName = results{1}(1).test;

    % To do: ensure all results have the same test type
    
    options = isp_interpretarguments(true, struct(...
        'evaluations', [], ...
        'xTickLabels', '', ...
        'modLabel', 'Modification', ...
        'distLabel', 'Distance measure', ...
        'linestyleoptions', struct, ...
        'integration', 'mean', ...
        'xTics', [], ...
        'xaxis', 'distancemeasure', ...
        'plotTitle', [upper(testName(1)) testName(2:end) ' results']), ...
                                     varargin{:});

    % Convert the structs to something tractable
    params=extractData(results, options.evaluations);


    %
    % Replace strings if necessary
    %
    if isfield(options, 'replace')
        r=options.replace;
        for fieldName={'iterationValue', 'modValues', 'dstName', 'evalCriterions'}
            f=fieldName{1};
            if ~isempty(params.(f))
                if iscell(params.(f))
                    for n=1:length(params.(f))
                        for m=1:size(r, 1)
                            if isequal(params.(f){n}, r{m,1})
                                params.(f){n} = r{m,2};
                                break
                            end
                        end
                    end
                else
                    for n=1:length(params.(f))
                        for m=1:size(r, 1)
                            if isequal(params.(f)(n), r{m,1})
                                params.(f)(n) = r{m,2};
                                break
                            end
                        end
                    end
                end
            end
        end
    end

    %
    % Set default values
    %
    
    if isempty(params.modValues)
        params.modValues = {''};
        titletext = @(ev, mo) sprintf('%s accuracy', ev);
    else
        titletext = @(ev, mo) sprintf('%s accuracy. %s: %s', ev, ...
                                      results{1}(1).modParameter, mo);
    end


    %
    % Reduce the amount of data to present
    %

    % Define the reduction
    switch results{1}(1).evalType
      case 'overlap'
        params.evalCriterions(~params.hasOverlapResult) = [];
        testRes = params.results_overlap(:,:,:,params.hasOverlapResult,:);
        testRes = flipdim(cumsum(flipdim(testRes, 5), 5), 5);
        testRes(:,:,:,:,1) = [];
        overlapAxVal{1} = sprintf('At least %d match', 1);
        for iVal = 2:size(testRes, 5)-1
            overlapAxVal{iVal} = sprintf('At least %d matches', iVal);
        end
        overlapAxVal{iVal+1} = sprintf('%d matching instr.', iVal+1);
      case {'nomix', 'reference'}
        testRes = params.results;
      case 'mixsets'
        testRes = params.results_all;
      otherwise
        error('Unknown evalType')
    end

    if size(testRes, 2) == 1
        %warning('Only one iteration - not plotting any statistics.')
        options.integration = 'all';
    end

    % Define how to reduce the amount of data
    useErrorBar = false;
    switch options.integration
      case 'mean'
        ite=@(x, y) mean(x, y);
        ax2 = {'Mean'};
      case 'stddev'
        if ~isfield(options.linestyleoptions, 'linestyle') && ...
                (~isfield(options, 'organization') || ...
                 strcmp(options.organization{2}, 'line l'))
            %ax2 = {'Mean', 'Mean $\pm$ std.dev', ''};
            ax2 = {'Mean', '', ''};
            ite=@(x, y) cat(y, mean(x, y), mean(x, y)+std(x, 0, y), ...
                        mean(x, y)-std(x, 0, y));
            %options.linestyleoptions.linestyle = {'-', '--', '--'};
            %options.linestyleoptions.linestyle = {'-', 'errorbar', 'errorbar'};
        else            
            ite=@(x, y) cat(y, mean(x, y), mean(x, y)+std(x, 0, y), ...
                        mean(x, y)-std(x, 0, y));
            ax2 = {'Mean', 'Mean + std.dev', 'Mean - std.dev'};
        end
        useErrorBar = true;

      case 'maxmin'
        ite=@(x, y) cat(y, max(x, [], y), min(x, [], y));
        ax2 = {'Maximum', 'Minimum'};
      case 'all'
        ite=@(x, y) x;
        ax2 = params.iterationValue;
      otherwise
        error('Invalid ''integration'' setting.')
    end

    % Do it
    testRes = ite(testRes, 2);


    %
    % Call the plot function
    %

    % Set common options
    ax1 = params.modValues;
    ax3 = params.dstName;
    ax4 = regexprep(params.evalCriterions, '$', ' accuracy');

    switch options.xaxis
      case 'modification'
        colors=linspace(0, 0.6, size(testRes, 1))'*[1 1 1];
      case 'distancemeasure'
        colors=linspace(0, 0.6, size(testRes, 3))'*[1 1 1];
      otherwise
        error('Invalid ''xaxis'' value.')
    end                

    if strcmp(results{1}(1).evalType, 'overlap')
        colors=linspace(0, 0.6, size(testRes, 5))'*[1 1 1];
    end


    linestyleoptions = struct('markeredgecolor', 'auto', ...
                              'markerfacecolor', 'none', ...
                              'color', colors);
% $$$     if strcmp(results{1}(1).evalType, 'mixsets')
% $$$         linestyleoptions.markeredgecolor = 'color';
% $$$         linestyleoptions.markerfacecolor = 'color';
% $$$     end

    options.linestyleoptions = isp_interpretarguments(linestyleoptions, ...
                                                      options.linestyleoptions);

    commonOptions = {'legendlocation', 'eastoutside', ...
                     'legendtype', 'type', ...
                     'ylim', [0 1], ...
                     'ylabel', 'Accuracy', ...
                     'figuretitle', options.plotTitle};

    % What the dimensions of testRes mean:
    %   1: modification (e.g. SNR)
    %   2: test set
    %   3: distance measure 
    %   4: evaluation criterion
    %   5: either modification axis 2 or 'overlap'

    moreOpts = struct;
    testResSize = [size(testRes) 1 1 1 1 1];
    switch results{1}(1).evalType
      case 'mixsets'
        switch options.xaxis
          case 'modification'
            organization={'xaxis', 'figure', 'line c|l', 'subploty', 'line m|e|f'};
          case 'distancemeasure'
            organization={'line c|l', 'figure', 'xaxis', 'subploty', 'line m|e|f'};
          otherwise
            error('Invalid ''xaxis'' value.')
        end                
        if useErrorBar, organization{2} = 'errorbar'; end
        if all(testResSize(find(strcmp(regexprep(organization,' .*',''),'line')))==1)
            moreOpts.legendlocation = 'off';
        end
        %isp_assert(strcmp(options.xaxis, 'distancemeasure'), ...
        %           'Invalid ''xaxis'' value.')
        plotargs = {testRes, commonOptions{:}, moreOpts, ...
              'organization', organization, ...
              'axisDescr', {options.modLabel, 'Test set', options.distLabel, ...
                            'Evaluation criterion', options.modLabel}, ...
              'axisValue', {ax1, ax2, ax3, ax4, ax1}, options};
        h = isp_plot(plotargs{:});

      case 'overlap'
        %   1: modification (e.g. SNR)
        %   2: test set
        %   3: distance measure 
        %   4: evaluation criterion
        %   5: overlap
        switch options.xaxis
          case 'modification'
            organization={'xaxis', 'line c|l', 'subploty', 'figure', 'line m'};
          case 'distancemeasure'
            organization={'figure', 'line c', 'xaxis', 'subploty', 'line l'};
          otherwise
            error('Invalid ''xaxis'' value.')
        end
        if useErrorBar, organization{2} = 'errorbar'; end
        linestyleoptions.color = linspace(0, 0.8, size(testRes, 2))'*[1 1 1];
        if all(testResSize(find(strcmp(regexprep(organization,' .*',''),'line')))==1)
            moreOpts.legendlocation = 'off';
        end
        plotargs = {testRes, commonOptions{:}, moreOpts, ...
              'organization', organization, ...
              'linestyleoptions', linestyleoptions, ...
              'axisDescr', {options.modLabel, 'Test set', options.distLabel, ...
                            'Evaluation criterion', 'Overlap'}, ...
              'axisValue', {ax1, ax2, ax3, ax4, overlapAxVal}, options};
        h = isp_plot(plotargs{:});

      case {'nomix', 'reference'}
        switch options.xaxis
          case 'modification'
            organization={'xaxis', 'line m', 'line cl', 'subploty'};
          case 'distancemeasure'
            organization={'line c|m', 'line l', 'xaxis', 'subploty'};
          otherwise
            error('Invalid ''xaxis'' value.')
        end
        if useErrorBar, organization{2} = 'errorbar'; end
        if all(testResSize(find(strcmp(regexprep(organization,' .*',''),'line')))==1)
            moreOpts.legendlocation = 'off';
        end
        plotargs = {testRes, commonOptions{:}, moreOpts, ...
              'organization', organization, ...
              'axisDescr', {options.modLabel, 'Test set', options.distLabel, ...
                            'Evaluation criterion'}, ...
              'axisValue', {ax1, ax2, ax3, ax4}, options};
        h = isp_plot(plotargs{:});
        
      otherwise
        error('fdsa')
    end


end % plotresultscomp_helper




function params = extractData(results, evalMask)
% Interprets the results in an cell array of structs of evaluation
% results.
%
% params is a struct with different information extracted from the results
% params.
%
% params.results(a,b,c,d) is the nearest neighbor accuracy of distance measure
% no. c wrt. evaluation criterion d (typically, d=1 is melody, d=2 is
% instrumentation) in test set b (e.g., sound font no.) with modification
% a. So,
%   a: modification (e.g. SNR)
%   b: test set
%   c: distance measure 
%   d: evaluation criterion
%
% params.results_same and params.results_diff have the same structure as params.res.
%    
% params.dstName: names of distance measures
% params.modValues: names of modifications measures
%
% params.overlap(a,b,c,d,e)
% e: nOverlaps

    if ~iscell(results)
        results={results};
    end

    if exist('evalMask', 'var') && ~isempty(evalMask)
        tests = results{1}(1).testsetname(evalMask);
    else
        tests = results{1}(1).testsetname;
    end

    nEvals = length(tests); % The number of ways to evaluate a result,

    overlap = reshape(nan, [1 1 1 1 1]);

    modValueCombinations = false;
    for iDstMsr=1:length(results)
        res=results{iDstMsr};
        params.dstName{iDstMsr} = res(1).distancename;
        if isfield(res(1), 'modValues')
            params.modValues = res(1).modValues;
        else
            params.modValues = [];
        end
        for jTest=1:nEvals
            for kIter = 1:numel(res)
                iterValues{kIter} = res(kIter).iterationValue;
                if any(strcmp(res(kIter).evalType, {'nomix', 'reference'}))
                    testRes(:,kIter, iDstMsr, jTest) = cat(1, res(kIter).(tests{jTest})(:).nnAccuracy);
                else
                    tmp = reshape(cat(1, res(kIter).(tests{jTest}).nnAccuracy), size(res(kIter).(tests{jTest})));
                    idx=size(tmp);
                    allRes(1:idx(1),kIter, iDstMsr, jTest, 1:idx(2)) = tmp;
                    sameRes(1:idx(1),kIter, iDstMsr, jTest) = diag(tmp);
                    if isscalar(tmp)
                        diffRes = zeros(0, kIter, iDstMsr, jTest);
                    else
                        diffRes(1:idx(1)*(idx(1)-1), kIter,iDstMsr,jTest) = ...
                            tmp(~eye(length(tmp)));
                    end
                    modValueCombinations = true;
                    testRes(1:prod(idx), kIter, iDstMsr, jTest) = cat(1, ...
                                          sameRes(:,kIter, iDstMsr, jTest), ...
                                          diffRes(:,kIter, iDstMsr, jTest));
                end
                if isempty(res(kIter).(tests{jTest})(1).overlapStatistics)
                    overlap(:,kIter, iDstMsr, jTest,:) = nan;
                else            
                    if isvector(res)
                        tmpMat = cat(1, res(kIter).(tests{jTest})(:).overlapStatistics.nearestNeighborAccuracy);
                        overlap(1:size(tmpMat,1),kIter, iDstMsr, ...
                                jTest,1:size(tmpMat,2)) = tmpMat;
                    else
                        error('Extracting results from the specified test configuration is not implemented.')
                    end
                end
            end
        end
    end


    if modValueCombinations && ~isempty(params.modValues)
        mv = params.modValues;
        idx=0;
        for n=1:length(mv)
            idx = idx + 1;
            params.modValues{idx} = mv{idx}; % somewhat unnecessary ... ;-)
        end
        for m=1:length(mv)
            for n=1:length(mv)
                if m==n, continue, end
                idx = idx + 1;
                params.modValues{idx} = [mv{m} ' vs. ' mv{n}];
            end
        end

    end

    params.evalCriterions = tests;
    params.results_overlap = overlap;
    params.iterationValue = iterValues;

    if exist('sameRes', 'var')
        params.results_same = sameRes;
        params.results_diff = diffRes;
        params.results_all = allRes;
    else
        params.results = testRes;
    end

    for n=1:size(overlap,4)
        hasOverlapResult(n) = any(reshape(~isnan(overlap(:,:,:,n,:)), [], 1));
    end
    params.hasOverlapResult = hasOverlapResult;
end % extractData

function str=any2str(arg)
    if iscell(arg)
        if numel(arg) > 1
            error('Cell structure with more than 1 element.')
        end
        arg=arg{1};
    end
    str=num2str(arg);
end
