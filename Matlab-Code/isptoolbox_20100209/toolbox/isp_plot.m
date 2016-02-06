%ISP_PLOT  Plot many-dimensional variable
%
% SYNTAX
%   h=isp_plot(data, options ...)
%   h=isp_plot(data, organization, options ...)
%
% DESCRIPTION
%   Plot 
%
% INPUT
%   data:
%     N-dimensional real matrix of data to plot.
%   organization:
%     Cell array containing N strings. The n'th element describes which
%     variable to change to illustrate variations along the n'th axis of
%     'data'. Possible values of each cell: 'xaxis', 'line', 'line cfg',
%     'subploty', 'subplotx'. In 'line cfg', "cfg" should be replaced by
%     a line change description.
%   options:
%     Structs with the following field names or  parameter, value pairs:
%     axisDescr:
%       Labels for the different dimensions of 'data'. The labels are
%       used to name the different axes (similar to xlabel).
%     axisIsOrdinal:
%       Array of booleans where the n'th element specifies whether the
%       labels of the n'th axis are ordinal. 
%     axisValue:
%       A cell array where the n'th element is an array with the values
%       of the different axes, as used in set(gca, 'XTickLabel', values).
%       So, axisValue{n} is an array with size(data, n) elements whose 
%       unit is described by axisDesc{n}.
%     figuretitle:
%       A function, cell array or string giving the figure
%       titles (not subplot titles). Default: The empty string.
%     legendincluded:
%       Array of booleans where element n specifies if the n'th the
%       legend should be displayed for the n'th line-dimension.
%     legendlocation:
%       'on', 'newfigure' or location.
%     legendorientation:
%       'horizontal' or 'vertical'.
%     legendtitle:
%       String, function or cell array.
%     legendbox:
%       Boolean specifying whether there should be a box around the
%       legend. Default: true.
%     legendtype:
%       'type' or 'all'.
%     linestyle:
%       A string configuring the line axes. It is recommended to leave
%       this option alone and specify line axes as part of 'organization'
%       instead.
%     linestyleoptions:
%       A struct configuring the line styles used. It is passed directly
%       to isp_linestyle.
%     plotparameters:
%       A string or cell array of strings giving additional arguments to
%       pass to the plot command. Default: Empty
%     subplotsep:
%       Separator between the description of the horizontal and vertical
%       subplots axes in the automatically generated subplot
%       titles. Default: ' - '
%     subplottitle:
%       A function, cell array or string  giving the subplot titles. If
%       unspecified, titles are automatically generated from 'axisValues'.
%     xlim, ylim:
%       X- and Y-axis limits. Default: [-inf inf], i.e., the plot
%       function's default.
%     xlog, ylog:
%       Booleans specifying whether to use logarithmic x- and y-axis,
%       respectively. Default: false.
%     ylabel:
%       The y-axis label. The x-axis label is taken from 'axisDesc'.
%
% OUTPUT
%   h:
%     Struct with an awful lot of handles to pretty much everything.
%
% EXAMPLE
%     data = rand(3, 4, 2, 5, 7);
%     organization = {'subplotx','line m','subploty','line l|c','xaxis'};
%     % A simple example
%     isp_plot(data, organization);
%     % Add some legends
%     descriptions = {'Age', 'Income', 'Sex', 'Height', 'Weight'};
%     axisValues = {[30 35 40], 10000:10000:40000, {'Man', 'Woman'}, ...
%                   150:10:190, 50:10:110};
%     isp_plot(data, organization, 'axisDescr', descriptions, ...
%         'axisValue', axisValues');
%     % Use customized subplot titles
%     isp_plot(data, organization, 'axisDescr', descriptions, ...
%              'axisValue', axisValues', 'subplottitle', ...
%              @(x,y) sprintf('%d year old %s', axisValues{1}(x), axisValues{3}{y}))
%
% SEE ALSO
%   isp_linestyle.
%
% HISTORY
%   Created January 2008 by Jesper H. Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [h]=isp_plot(data, varargin)

    nDim = ndims(data);

    validDimensions = {'xaxis', 'line', 'subploty', 'subplotx'};

    iVararg=1;
    if length(varargin)>=1 && iscell(varargin{1})
        options.organization = varargin{1};
        iVararg = iVararg + 1;
    else
        options.organization = {'xaxis', 'line cl', 'subploty', 'subplotx'};
    end

    for n=1:nDim
        options.axisDescr{n} = ['Axis ' num2str(n)];
        options.axisValue{n} = 1:size(data, n);
    end

    options.ylabel = '';
    options.linestyle = '';
    options.linestyleoptions = struct();
    options.plotparameters = '';
    options.xlog = false;
    options.ylog = false;
    options.xlim = [-inf inf];
    options.ylim = [-inf inf];
    options.figuretitle = nan;
    options.subplottitle = nan;
    options.legendbox = true;
    options.legendorientation = 'vertical';
    options.legendtitle = nan;
    options.legendlocation = 'newfigure';
    options.legendtype = 'type';
    options.legendincluded = [];
    options.subplotsep = ' - ';
    options.axisIsOrdinal(1:nDim) = true;
    options = isp_interpretarguments(options, varargin{iVararg:end});

    if length(options.organization) < nDim
        options.organization(end+1:nDim) = {'figure'};
    end
    options.organization(end+1:nDim) = [];

    xaxDim = nDim+1;
    lineDim = nDim+1;
    spxDim = nDim+1;
    spyDim = nDim+1;
    plotDim = nDim+1;
    errorBarDim = nDim+1;
    
    for n=1:nDim
        plotSpec = options.organization{n};
        if strcmp(plotSpec(1:4), 'line') && ...
                length(plotSpec) > 4
            if n==1, options.linestyle=''; end
            options.linestyle = [options.linestyle ' ' plotSpec(5:end)];
            plotSpec = 'line';
        end

        switch plotSpec
          case 'xaxis'
            xaxDim = n;
          case 'line'
            if isequal(lineDim, nDim+1), lineDim = n;
            else lineDim(end+1) = n; end
          case 'subplotx'
            spxDim = n;
          case 'subploty'
            spyDim = n;
          case 'errorbar'
            errorBarDim = n;
          case 'figure'
            if isequal(plotDim, nDim+1), plotDim = n;
            else plotDim(end+1) = n; end
          otherwise
            error(['Invalid dimension specification ''' plotSpec '''.'])
        end
    end

    options.legendincluded = isp_nonempty(options.legendincluded, ...
                                          true(size(lineDim)));
    s=[size(data) 1];
    nX = s(xaxDim);
    nSpy = s(spyDim);
    nSpx = s(spxDim);
    nLines = prod(s(lineDim));
    nPlots = prod(s(plotDim));
    nErrorBars = s(errorBarDim);
    
    if nErrorBars > 3
        error('Data cannot have size larger than 3 in the errorBar dimension.')
    end
    
    if isequalwithequalnans(nan, options.figuretitle)
        if length(plotDim) == 1 && plotDim > nDim
            options.figuretitle = '';
        else
            options.figuretitle = @(varargin) '';
        end
    end
    
    if isequalwithequalnans(nan, options.subplottitle)
        if spxDim > nDim
            xspText= @(x,y) '';
        else xspText = @(kSpx, jSpy) any2str(options.axisValue{spxDim}(kSpx));
        end
        if spyDim > nDim
            yspText= @(x,y) '';
        else yspText = @(kSpx, jSpy) any2str(options.axisValue{spyDim}(jSpy));
        end
        % Sigh, I miss (exp) ? a:b ...
        if spyDim <= nDim && spxDim <= nDim
            tmp = options.subplotsep;
        else tmp = '';
        end
        options.subplottitle = @(kSpx, jSpy) ...
            [xspText(kSpx, jSpy) tmp yspText(kSpx, jSpy)];
    end

    if isequalwithequalnans(nan, options.legendtitle)
        if length(lineDim) == 1 && lineDim > nDim
            options.legendtitle = '';
            options.legendlocation = 'off';
        else
            for n=1:length(lineDim)
                if numel(options.axisValue{lineDim(n)}) < s(lineDim(n))
                    error(['Fewer elements in ''axisValue'' than in ' ...
                           '''data'' in dimension ' num2str(lineDim(n)) '.'])
                end
            end
            switch options.legendtype
              case 'all'
                tmpfunc = @(varargin) '';
                for n=1:length(lineDim)
                    f = @(x) [options.axisDescr{lineDim(n)} ' ' ...
                              any2str(options.axisValue{lineDim(n)}(x))];
                    tmpfunc = @(varargin) [tmpfunc(varargin{:}) f(varargin{n})];
                end
                options.legendtitle = tmpfunc;
              case 'type'
                A=options.axisValue(lineDim);
                options.legendtitle = @(varargin) ...
                   any2str(A{cat(1, varargin{:})>0}(max(cat(1, varargin{:}))));
              otherwise                
                error('Invalid value of legendtype.')
            end
        end
    end

    options.linestyle = isp_nonempty(options.linestyle, 'l');

    for iPlot = 1:nPlots
        plotIdx={};
        [plotIdx{1:length(plotDim)}] = ind2sub(s(plotDim), iPlot);
        h.figure(plotIdx{:}) = figure;
        plotIdxMat = cat(2, plotIdx{:});
        for jSpy = 1:nSpy

            for kSpx = 1:nSpx
                % When resizing subplots with external legends, Matlab
                % handles resizing incorrectly (at least in
                % 2007b). Therefore we only create a subplot if
                % necessary.
                if nSpy~=1 || nSpx~=1
                    h.subplot(jSpy,kSpx,iPlot) = subplot(nSpy, nSpx, ...
                                                         kSpx+(jSpy-1)*nSpx);
                end        
                
                for ting=1:2
                    legText = [];
                    for lLine = 1:nLines
                        lineIdx={};
                        [lineIdx{1:length(lineDim)}] = ind2sub(s(lineDim), lLine);
                        lineIdxMat = cat(2, lineIdx{:});
                    
                        idx=ones(nX, nDim+1);
                        idx(:, spyDim) = jSpy;
                        idx(:, spxDim) = kSpx;
                        idx(:, lineDim) = repmat(lineIdxMat, nX, 1);
                        idx(:, plotDim) = repmat(plotIdxMat, nX, 1);
                        idx(:, xaxDim) = (1:nX)';
                        idx(:, errorBarDim) = 1;
                        idc = num2cell(idx, 1);
                        val = data(sub2ind(s, idc{:}));
                        if length(options.plotparameters) > 1
                            par = options.plotparameters{iLine};
                        else
                            par = options.plotparameters;
                        end
                        if ischar(par)
                            par = {par};
                        end
                        
                        linestyle{lLine} = ...
                            isp_linestyle(lineIdxMat, options.linestyle, options.linestyleoptions);
                        par(end+1:end+numel(linestyle{lLine})) = ...
                            linestyle{lLine};
                        
                        myhold
                        if ting==2
                            if iscell(options.axisValue{xaxDim}) || ...
                                    options.axisIsOrdinal(xaxDim)
                                xval = 1:length(options.axisValue{xaxDim});
                                plot(xval, val, ...
                                     par{:});
                                set(gca, 'XTick', 1:length(options.axisValue{xaxDim}))
                                set(gca, 'XTickLabel', options.axisValue{xaxDim})
                            else
                                xval = options.axisValue{xaxDim};
                                plot(xval, val, par{:});
                            end

                            if nErrorBars > 1
                                idx(:, errorBarDim) = 2;
                                idc = num2cell(idx, 1);
                                errValLow = data(sub2ind(s, idc{:}));

                                switch nErrorBars
                                  case 2
                                    errValUpp = val + errValLow;
                                    errValLow = val - errValLow;
                                  case 3
                                    idx(:, errorBarDim) = 3;
                                    idc = num2cell(idx, 1);
                                    errValUpp = data(sub2ind(s, idc{:}));
                                end  

                                errWidth = 0.1*mean(abs(diff(xval)));
                                for n=1:length(xval)
                                    for newVal = {errValUpp, errValLow}
                                    plot(xval([n n]), [val(n) newVal{1}(n)], ...
                                         par{:}, 'Marker', 'none');
                                    plot(xval(n)+[-1 1]*errWidth, ...
                                         newVal{1}([n n]), par{:}, ...
                                         'Marker', 'none', 'LineStyle', '-');
                                    end
                                end
                            end


                            if options.xlog
                                set(gca, 'XScale', 'log')
                            end
                            if options.ylog
                                set(gca, 'YScale', 'log')
                            end
                        end

                        if strcmp(options.legendtype, 'all')
                            if isa(options.legendtitle, 'function_handle')
                                legText{lLine} = options.legendtitle(lineIdx{:});
                            elseif iscell(options.legendtitle)
                                legText{lLine} = options.legendtitle{lineIdx{:}};
                            else legText{lLine} = options.legendtitle;
                            end
                        end
                    end

                    if ting==1
                        switch options.legendlocation
                          case {'on', true}
                            [h.legend{jSpy,kSpx,iPlot}{1:4}] = mylegend(...
                                options, s(lineDim), legText, ...
                                'Orientation', options.legendorientation);

                          case 'off'
                          case 'newfigure'
                          otherwise
                            h.legend{jSpy,kSpx,iPlot} = cell(1,4);
                            [h.legend{jSpy,kSpx,iPlot}{1:4}] = mylegend(...
                                options, s(lineDim), legText, ...
                                'Location', options.legendlocation, ...
                                'Orientation', options.legendorientation);
                        end
                        
                    end
                end



                xlim(options.xlim)
                ylim(options.ylim)
                switch class(options.subplottitle)
                  case 'function_handle'
                    title(options.subplottitle(kSpx, jSpy))
                  case 'cell'
                    title(options.subplottitle{kSpx, jSpy})
                  otherwise
                    % Let's hope it's a string ...
                    title(options.subplottitle)
                end
                xlabel(options.axisDescr{xaxDim})
                ylabel(options.ylabel)
            end
        end
        switch class(options.figuretitle)
          case 'function_handle'
            isp_subplotlabel('title', options.figuretitle(plotIdx{:}));
          case 'cell'
            isp_subplotlabel('title', options.figuretitle{plotIdx{:}});
          otherwise
            % Let's hope it's a string ...
            if ~isempty(options.figuretitle)
                isp_subplotlabel('title', options.figuretitle);
            end
        end
        if nSpy > 1
            h.subplotlabely{iPlot} = ...
                isp_subplotlabel('y', options.axisDescr{spyDim});
        end
        if nSpx > 1
            h.subplotlabelx{iPlot} = ...
                isp_subplotlabel('x', options.axisDescr{spxDim});
        end
    end

    switch options.legendlocation
      case 'on'
      case 'off'
      case 'newfigure'
        figure
% $$$         for n=1:nLines
% $$$             plot(1, nan, linestyle{n}{:})
% $$$             myhold
% $$$         end
        [h.legend{1}{1:4}] = mylegend(options, s(lineDim), legText, ...
                                      'Location', 'North', 'Orientation', ...
                                      options.legendorientation);
        axis off
      otherwise
        % 'Location' value for legend
    end

end



function varargout=mylegend(options, sizeLineDim, legText, varargin)

    % First handle the simple case ...
    if strcmp(options.legendtype, 'all')
        [varargout{1:4}] = legend(legText{:}, varargin{:});
        if options.legendbox
            legend boxon
        else
            legend boxoff
        end
        return
    end

    if ~strcmp(options.legendtype, 'type')
        error('Invalid value of legendtype.')
    end

    
    linopts.color = [0 0 0];
    linopts.linestyle = {'-'};
    linopts.markerstyle = {'.'};
    if isfield(options.linestyleoptions, 'linewidth') && ...
            length(options.linestyleoptions.linewidth) == 1
        linopts.linewidth = options.linestyleoptions.linewidth;
    else
        linopts.linewidth = .5;
    end
    if isfield(options.linestyleoptions, 'markersize') && ...
            length(options.linestyleoptions.markersize) == 1
        linopts.markersize = options.linestyleoptions.markersize;
    else
        linopts.markersize = 6;
    end
    linopts.markeredgecolor = {'k'};
    linopts.markerfacecolor = {'none'};

    legText = {};

    axesWithLegendText = find(options.legendincluded);
    for iAxis = 1:length(axesWithLegendText)
        kDim = axesWithLegendText(iAxis);
        idx=ones(size(sizeLineDim));
        idx2=zeros(size(sizeLineDim));
        nonEmpty = true(sizeLineDim(kDim),1);
        for jDim = 1:sizeLineDim(kDim);
            idx(kDim) = jDim;
            idx2(kDim) = jDim;
            lineIdx = num2cell(idx);
            lineIdx2 = num2cell(idx2);
            if isa(options.legendtitle, 'function_handle')
                legText{iAxis}{jDim,1} = options.legendtitle(lineIdx2{:});
            elseif iscell(options.legendtitle)
                legText{iAxis}{jDim,1} = options.legendtitle{lineIdx{:}};
            else
                legText{iAxis}{jDim,1} = options.legendtitle;
            end
            nonEmpty(jDim) = ~isempty(legText{iAxis}{jDim,1});
            [linedef{iAxis,1}{jDim,1}, dummy, change{iAxis}(jDim,1:2)] = ...
                isp_linestyle(idx, options.linestyle, ...
                                    options.linestyleoptions, 'initialvalues', ...
                                    linopts, 'selecteddimensions', kDim);
        end
        legText{iAxis}{jDim+1,1} = '';
        linedef{iAxis}{jDim+1,1} = {};
        change{iAxis}(jDim+1,1:2) = [false false];
        legText{iAxis}(~nonEmpty) = [];
        linedef{iAxis}(~nonEmpty) = [];
        change{iAxis}(~nonEmpty,:) = [];
    end
    ld=linedef;
    lt=legText;
    linedef = cat(1, linedef{:});
    legText = cat(1, legText{:});
    change = cat(1, change{:});








    %nEntries = sum(sizeLineDim(axesWithLegendText)) + numel(axesWithLegendText);
    nEntries = length(linedef)-1;

    % Dummy plots are needed if the legend is in a separate
    % figure. It shouldn't hurt if we aren't.
    hold all
    
    for n=1:nEntries
        % I much prefered using nan, but unfortunately that was
        % unreliable

        if change(n,1) % Line
            args=linedef{n}(1:6);
        else
            args={'linestyle', 'none'};
        end
        if change(n,2) % marker
            args(end+1:end+length(linedef{n})-6)=linedef{n}(7:end);
        end
        plot(mean(xlim), inf, args{:});
    end
    [varargout{1:4}] = legend(legText{1:end-1}, varargin{:});
    objh = varargout{2};
    if options.legendbox
        legend boxon
    else
        legend boxoff
    end
end



function str=any2str(arg)
    if iscell(arg)
        if numel(arg) > 1
            error('Cell structure with more than 1 element.')
        end
        arg=arg{1};
    end
    str=num2str(arg);
end

function myhold
    if exist('OCTAVE_VERSION')
        hold on
    else
        hold all
    end
end