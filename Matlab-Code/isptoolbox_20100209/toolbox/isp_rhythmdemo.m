%% ISP_RHYTHMDEMO  Demonstration of isp_tirhythm and the evaluation framework
% Test different variations of the tempo-insensitive rhythm feature

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


datapath = 'evaluatedata';
commonOptions = struct('experiment', 'ballroom', ...
                       'ismir2004ballroompath', './', ...
                       'dataPath', datapath);

%%
% We use the cacheRawData feature of the isp_tirhythm feature. This makes it
% possible to do the computationally heavy tempo extraction once and for
% all, and postpone the lightweight aggregation into the compact
% tempo-insensitive feature until the distance computation stage. This
% allows us to reuse the tempo extraction.

%%
% Define a dummy distance measure that only extracts the tempo
dummyRhythMeasure = isp_tirhythm;
dummyRhythMeasure.options.cacheRawData = true;

%%
% Extract the features
isp_evaluate(dummyRhythMeasure, commonOptions, 'distribute', 1);
isp_runjob(fullfile(datapath, 'distributedfeatures'))

%%
% Define the real features that we actually want to test
distanceMeasures = {};
distanceMeasures{end+1} = isp_tirhythm;
% The paper where the tempo-insensitive rhythmic patterns are introduced ("A
% TEMPO-INSENSITIVE REPRESENTATION OF RHYTHMIC PATTERNS" by J. H. Jensen,
% M. G. Christensen and S. H. Jensen), uses 60 bands. The current default is
% 45 bands, since accuracy seems nearly identical at a lower
% dimensionality. Uncomment the next line to duplicate the results from the
% paper.
%distanceMeasures{end}.options.nBands = 60;
distanceMeasures{end}.name = 'Tempo insensitive';

% The non-tempoinsensitive version
distanceMeasures{end+1} = isp_tirhythm;
distanceMeasures{end}.name = 'Linear';
distanceMeasures{end}.options.variation='raw';

% A re-implementation of the feature in ISMIR 2007 paper "From Rhythm
% Patterns to Perceived Tempo" by K. Seyerlehner, G. Widmer and D. Schnitzer.
distanceMeasures{end+1} = isp_tirhythm;
distanceMeasures{end}.name = 'Seyerlehner'
distanceMeasures{end}.options.variation='rawseyerlehner';
distanceMeasures{end}.options.tMin=0;
distanceMeasures{end}.options.offset=0;

%%
% Compute distances and evaluate the results
for n=1:length(distanceMeasures)
    dstMsr = distanceMeasures{n};
    dstMsr.options.cacheRawData = true;
    results{n} = isp_evaluate(dstMsr, commonOptions, 'distribute', 3);
end

%%
% Plot and print the results

isp_plotresults(cat(1, results{:}));

fprintf('\n\n');
for n=1:length(distanceMeasures)
    fprintf('%s: Accuracy %f. Accuracy when discarding songs of similar tempo: %f\n', results{n}{1}.distancename, results{n}{1}.split.nnAccuracy, results{n}{1}.difftemposplit.nnAccuracy);
end
fprintf('\n\n');


%%
% Just to show how to extract the compact features directly, we plot the raw
% features for BallroomData/ChaChaCha/Albums-Cafe_Paradiso-05.wav for the
% three distance measures we have defined. In this case, we do not want
% any intermediate feature, so we need to set 'cacheRawData' to false

file = fullfile(isp_toolboxpath, 'Loveshadow - The_Acorns. Seedin Time in The Oak Room - excerpt.mp3');
figure
for n=1:length(distanceMeasures)
    dstMsr = distanceMeasures{n};
    dstMsr.options.cacheRawData = false;
    feature = isp_extractfeature(file, dstMsr);
    subplot(length(distanceMeasures), 1, n)
    bar(feature)
    title(dstMsr.name)
end
