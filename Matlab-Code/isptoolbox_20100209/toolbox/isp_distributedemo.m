%% isp_distributedemo  Demonstration of isp_distribute.
% In the following, we demonstrate how to distribute stuff to multiple
% computers. We have a script named isp_distributedemo_helper1, which
% calls the (slow) function  isp_distributedemo_helper2.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.

type isp_distributedemo_helper1

%%
%
type isp_distributedemo_helper2

%%
% During the first iteration, we distribute jobs
isp_distribute('start')
isp_distributedemo_helper1

%%
% We see that neither res nor callNo exist yet, and .
whos res callNo

%%
% Second time, we run the jobs. For this demonstration, we
% "cheat" and just use MATLAB's parfor to run 2 parallel threads,
% but it could just as well have been on separate computers.
% When the functions we execute are slow, it is almost inevitable that
% we will get some error messages when two threads attempt to access the
% same job file simultaneously.
matlabpool open 2
parfor (m=1:4)
    isp_runjob
end
matlabpool close

%%
% Third time we want to read the results
isp_distribute('read')
isp_distributedemo_helper1

res
callNo
