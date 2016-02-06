%% ISP_EVALUATEDEMO2  Demonstration of how to use the evaluation framework
% In the following, we evaluate how well the isp_tichroma distance
% measure behaves when songs are identical except for tempo.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


dstMsr = isp_tichroma;
datapath = './evaluatedata';

%%
% We only run the test with four different instruments and three
% melodies. For any serious use, use the default values of 30 instruments
% and 30 melodies.
tic
res = isp_evaluate(dstMsr, ...
             'distribute', 'ask', ...
             'experiment', {'instrumentmelody', 'duration'}, ...
             'dataPath', datapath, ...
             'nInstruments', 4, ...
             'nMidifiles', 3);
toc

%%
% Finally, we plot the results. In Figure 1, we see that the melody
% recognition accuracy is 100%, and the instrument recognition accuracy is
% 0, which is what we would except for a distance measure developed for
% cover song recognition. In Figure 2, we see that when the length of a song
% is 3/4 of the original, the melody recognition accuracy drops to somewhere
% between 60% and 70%. Again, with only 4 instruments and 3 melodies in the
% test set, these results should be taken with a grain of salt.
isp_plotresults(res, 'xaxis', 'modification');
