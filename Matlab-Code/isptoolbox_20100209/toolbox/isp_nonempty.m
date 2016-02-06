%ISP_NONEMPTY  Return the first non-empty input argument
%
% SYNTAX
%   arg=isp_nonempty(input1, input2, ...)
%
% DESCRIPTION
%   Returns the first non-empty input argument. It is intended for use in
%   situations where an empty value should be replaced by a default
%   value. It is particularly useful with isp_interpretarguments when
%   some default values depend on other input arguments.
%
% INPUT
%   input1, input2, ...:
%     Any inputs.
%
% OUTPUT
%   arg:
%     The first non-empty input argument. If all input arguments are
%     empty, the last argument is returned.
%
% EXAMPLE
%     function myplot(x, y, varargin)
%       % Allow user to specify empty x:
%       x = isp_nonempty(x, 1:length(y));
%       plot(x,y)
%
% SEE ALSO
%   isp_interpretarguments.
%
% HISTORY
%   Created by Jesper H. Jensen January 2008.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function arg=isp_nonempty(varargin)
    for n=1:nargin-1
        if ~isempty(varargin{n})
            arg=varargin{n};
            return
        end
    end
    arg=varargin{end};
end
