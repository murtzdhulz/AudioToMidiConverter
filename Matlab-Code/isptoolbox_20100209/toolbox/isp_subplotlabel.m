%ISP_SUBPLOTLABEL  Plot per-figure labels 

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.

function [ax,h]=isp_subplotlabel(varargin)
    if nargin <= 1 || isempty(varargin{1})
        return
    end


    options = isp_interpretarguments(struct('x', '', 'y', '', 'title', ''), ...
                                     varargin{:});

    ax=struct;
    h=struct;

    if ~exist('suplabel')
        warning(['Cannot find the function suplabel. You need to download ' ...
                 'this separately.'])
    else

        if ~isempty(options.x)
            [ax.x, h.x] = suplabel(options.x, 'x');
        end

        if ~isempty(options.y)
            [ax.y, h.y] = suplabel(options.y, 'y');
        end

        if ~isempty(options.title)
            [ax.title, h.title] = suplabel(options.title, 't');
        end
    end
end
