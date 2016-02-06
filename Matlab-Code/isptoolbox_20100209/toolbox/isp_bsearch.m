%ISP_BSEARCH  Binary search
%
% SYNTAX
%   idx=isp_bsearch(borders, val)
%
% DESCRIPTION
%   Find the index 'idx' in 'borders' such that
%     borders(idx) <= val < borders(idx+1)
%
% INPUT
%   borders:
%     Sorted array of border values.
%   val:
%     The value to search for. If 'val' is an array, 'idx' will also be
%     an array of the same size.
%
% OUTPUT
%   idx:
%     Index of the element or interval in 'borders' that match 'val'. If
%     val < borders(1), idx=0 is returned, and if borders(end) <= val, 
%     idx=length(borders)+1 is returned.
%
% EXAMPLE
%   borders=round(sort(rand(1,10)*100))
%   isp_bsearch(borders, 0.6)
%
% HISTORY
%   Created by Jesper H. Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function idx=isp_bsearch(borders, val)

    idx=zeros(size(val));

    for n=1:numel(idx)
        maxIndex=length(borders);
        minIndex=1;
        while maxIndex-minIndex ~= 1
            middleIndex=ceil(0.5*(maxIndex + minIndex));
            if val(n) >= borders(middleIndex)
                minIndex=middleIndex;
            else
                maxIndex=middleIndex;
            end
        end
        idx(n)=maxIndex-1;
    end

    idx(val >= borders(end)) = maxIndex+1;
    idx(val < borders(1)) = 0;
end