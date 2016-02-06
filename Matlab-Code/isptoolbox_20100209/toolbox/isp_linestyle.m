%ISP_LINESTYLE  Create options for line and marker styles.
%
% SYNTAX
%   [param, options] = isp_linestyle(idx, varargin)
%   [param, options] = isp_linestyle(idx, order, varargin)
%
% DESCRIPTION
%   This function makes it easy to create custom linestyles when MATLAB's
%   axes' ColorOrder and LineStyleOrder properties are not flexible enough.
%   It can be used to just create a lot of different linestyles, or it can be
%   used when e.g. the lines' color should depend on one parameter, the style
%   on another, and the marker type on a third parameter. Although isp_linestyle
%   can be called directly, in many cases it is easier to use isp_plot.
%
% INPUT
%   idx:
%     Array of length N, where N is the number of parameters the lines'
%     appearance shall depend on (see the example section). The n'th element
%     specifies 
%   order:
%     String configuring the line styles. It should contain numel(idx)
%     space-delimited 'words' that configures each dimension. Each 'word'
%     consists of one or more of the following letters:
%     'c':
%       Color
%     'l':
%       linestyle
%     'm':
%       marker
%     'w':
%       linewidth
%     's':
%       markersize
%     'e':
%       markeredge
%     'f':
%       markerface
%     If the letters of a 'word' are delimited by '|', the properties will
%     change simultaneously, e.g. 'c|l' will change both color and line
%     style simultaneously. If there is no delimiter between letters, the
%     properties will be changed in turn, e.g. 'cl' will only change line
%     style after having cycled through all color combinations, after which
%     it once again cycles through all color combinations. As an example,
%     if idx is 3-dimensional, 'order' could be e.g. 'cl m e|f', which
%     specifies that for the first dimension, we cycle through all colors
%     for each change of line style, similar to what MATLAB usually does
%     for plots. For the second dimension, we cycle through different marker
%     styles, and finally we cycle through the marker colors, this time
%     changing the edge and face color simultaneously.
%     For 1 to 4 dimensional idx, the default values are 'c', 'c l', 
%     'c l m', and 'c l m w', respectively.
%   options ...:
%     Structs or field/value pairs with the following field names:
%       color:
%         An Nx3 matrix with colors to switch between.
%         Default: get(0, 'DefaultAxesColorOrder')
%       initialvalues:
%         Struct with linestyle, linewidth, markerstyle, markersize,
%         markeredgecolor and markerfacecolor fields, respectively, that
%         have the same form as the options with the same names. The first
%         element of the 'inittialvalues' fields are used when nothing is
%         specified by 'order', possibly as a consequence of the values of
%         'selecteddimensions'. Neither 'initialvalues' nor
%         'selecteddimensions' should probably be used directly. They are
%         used for rather specialized interaction with isp_plot.
%       linestyle:
%         Cell array with linestyles to switch between.
%         Default: {'-', '--', ':', '-.'}
%       linewidth:
%         Line widths to switch between. Default: [0.5 2 5];
%       markeredgecolor:
%         Cell array with marker edge colors to switch between. 
%         Default: {'auto'};
%       markerfacecolor:
%         Cell array with marker face colors to switch between.
%         Default: {'none'};
%       markersize:
%         Marker sizes to switch between. Default: 6 (i.e., the same is
%         always used)
%       markerstyle:
%         Cell array with linestyles to switch between.
%         Default: {'.', 'o', '+', '*', 's', '^'}
%       order:
%         Alternative way to specify the order string.
%       suppresswarning:
%         Boolean telling whether to issue a warning if there are too few
%         combinations of line style, color etc. to give unique combinations.
%       selecteddimensions:
%         Only return the plotting arguments that are due to the dimensions 
%         given by this argument. Default: 1:numel(idx), i.e., return all
%         arguments.
%     
% OUTPUT
%   param:
%     The parameters to pass to the plot routine. See how to use them in 
%     the example section.
%   options:
%     Struct with input options supplemented with default values.
%
% EXAMPLE
%   Let's say we have some simulation results where we have changed
%   parameter x, y and z between 3, 2 and 4 different values each,
%   and that the results are in the [nx ny nz]-dimensional variable 'results'.
%   We plot the z values along the x-axis, but now we want changes in x and y
%   to be reflected by the lines' appearances. Below, we let color and line
%   style reflect the x axis, and the marker style reflect the y axis.
%
%     % Create artificial data
%     nx=3; ny=2; nz=4;
%     results = rand([nx ny nz]);
%     % Do the actual plotting.
%     clf; hold all
%     for ix = 1:nx
%       for iy = 1:ny
%         plotArgs = isp_linestyle([ix iy], 'l|c m');
%         plot(squeeze(results(ix, iy, :)), plotArgs{:})
%       end
%     end
%
% SEE ALSO
%   isp_plot, plot.
%
% HISTORY
%   2007: Created by Jesper H. Jensen 

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [param, options, change] = isp_linestyle(idx, varargin)
    nDims = numel(idx);

    options.color = get(0, 'DefaultAxesColorOrder');
    %options.linestyle = {'none', '-', '--', ':', '-.'};
    options.linestyle = {'-', '--', ':', '-.'};
    %options.markerstyle = {'none', '.', 'o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h'};
    options.markerstyle = {'.', 'o', '+', '*', 's', '^'};
    options.markersize = 6;
    options.linewidth = [0.5 2 5];
    options.markeredgecolor = {'auto'};
    options.markerfacecolor = {'none'};
    options.suppresswarning = false;
    options.initialvalues = [];    
    options.selecteddimensions = [];

    switch(length(idx))
      case 1
        options.order = 'c';
      case 2
        options.order = 'c l';
      case 3
        options.order = 'c l m';
      case 4
        options.order = 'c l m w';
      otherwise
        options.order = '';
    end

    if nargin >= 1 && ischar(varargin{1}) && ~isfield(options, varargin{1})
        options.order = varargin{1};
        iVararg = 2;
    else
        iVararg = 1;
    end

    options = isp_interpretarguments(options, varargin{iVararg:end});
   

    if isequal(lower(options.markeredgecolor), 'color')
        options.markeredgecolor = options.color;
    end

    if isequal(lower(options.markerfacecolor), 'color')
        options.markerfacecolor = options.color;
    end

    if ~iscell(options.markerfacecolor)
        options.markerfacecolor = num2cell(options.markerfacecolor, 2);
    end

    if ~iscell(options.markeredgecolor)
        options.markeredgecolor = num2cell(options.markeredgecolor, 2);
    end

    tmp = isp_nonempty(options.initialvalues, options);
    color = tmp.color(1,:);
    linestyle = tmp.linestyle{1};
    linewidth = tmp.linewidth(1);
    marker = tmp.markerstyle{1};
    markersize = tmp.markersize(1);
    markeredge = tmp.markeredgecolor{1};
    markerface = tmp.markerfacecolor{1};
    
    layout=split(options.order);

    isp_assert(length(layout) == nDims, ...
               'Invalid number of dimensions in ''order'' string.');

    desc='clwmsef';
    nExpressions(desc) = [size(options.color, 1)
                          length(options.linestyle)
                          length(options.linewidth)
                          length(options.markerstyle)
                          length(options.markersize)
                          size(options.markeredgecolor, 1)
                          size(options.markerfacecolor, 1)];

    tmp = isp_nonempty(options.selecteddimensions, 1:nDims);

    lineChange = false;
    markerChange = false;

    for n=tmp
        combis = split(layout{n}, '|');
        overflowWarning = ~options.suppresswarning;
        for m=1:length(combis)
            b = nExpressions(combis{m});
            c=cumprod(b);
            overflowWarning = overflowWarning & (idx(n) > c(end));
            indices = 1+mod([idx(n)-1 floor((idx(n)-1)./c(1:end-1))], b);
            for k=1:length(indices)
                lineChange = lineChange | any(combis{m}(k) == 'clw');
                markerChange = markerChange | any(combis{m}(k) == 'msef');
                switch(combis{m}(k))
                  case 'c', color = options.color(indices(k),:);
                  case 'l', linestyle = options.linestyle{indices(k)};
                  case 'm', marker = options.markerstyle{indices(k)};
                  case 'w', linewidth = options.linewidth(indices(k));
                  case 's', markersize = options.markersize(indices(k));
                  case 'e', markeredge = options.markeredgecolor{indices(k)};
                  case 'f', markerface = options.markerfacecolor{indices(k)};
                  otherwise, error(['Invalid format character ' c{n} '.'])
                end
            end % k
        end % m
        if overflowWarning
            warning(['Too many values in dimension ' num2str(n) ' to plot them all differently'])
        end
    end % n
            
    param = {'Color', color, ...
             'LineStyle', linestyle, ...
             'LineWidth', linewidth, ...
             'Marker', marker, ...
             'MarkerSize', markersize, ...
             'MarkerEdgeColor', markeredge, ...
             'MarkerFaceColor', markerface};
    change = [lineChange markerChange];
end

function c=split(str, varargin)
    iTok=1;
    c={};
    while ~isempty(str)
        [t,str] = strtok(str,varargin{:});
        c{iTok} = t;
        iTok = iTok + 1;
    end
end
