%ISP_DISTRIBUTEDEMO_HELPER1  Helper function for isp_distributedemo

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


for n=1:10
    isp_distribute('[res(n), callNo(n)] = isp_distributedemo_helper2(n)');
end

