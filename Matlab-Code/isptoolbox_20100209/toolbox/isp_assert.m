%ISP_ASSERT  Print error message if assertion does not hold.
%
% SYNTAX
%   isp_assert(expression, errorMessage)
%
% DESCRIPTION
%   Throw an error if expression is not true.
%
% INPUT
%   expression:
%     An expression that evaluates to True if everything is working as
%     expected.
%   errorMessage:
%     Message to print if expression evaluates to False.
%
% EXAMPLE
%     isp_assert(mod(N, 2)==0, 'Even N expected')
%
% HISTORY
%   Created by Jesper H. Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function isp_assert(varargin)

    if nargin < 1
        error('At least one input required.')
    end

    switch class(varargin{1})
      case 'logical'
        if ~varargin{1}
            if nargin>1
                error(varargin{2:end})
            else
                error('Assertion does not hold.')
            end
        end
      otherwise
        error('Invalid input argument.')
    end

end