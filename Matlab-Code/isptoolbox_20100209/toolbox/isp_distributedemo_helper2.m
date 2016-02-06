%ISP_DISTRIBUTEDEMO_HELPER2  Helper function for isp_distributedemo

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


% This function simulates a very slow function. In reality, it does
% nothing but to return the input variable and the number of times it has
% been called.
function [res, callNoOut] = isp_distributedemo_helper2(n)

    persistent callNo
    if isempty(callNo), callNo = 0; end
    callNo = callNo + 1;
    callNoOut = callNo;

    % In real life, this would have been a very slow function
    res = n;
    pause(2);

end
