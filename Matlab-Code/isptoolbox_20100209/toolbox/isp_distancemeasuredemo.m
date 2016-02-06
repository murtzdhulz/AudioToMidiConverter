%% ISP_DISTANCEMEASUREDEMO  Demonstrate various aspects of music distance measures
% In the following, we create a music distance measure and extracts
% features from a song, compute distances and convert it to/from an array
% of doubles
% 

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


%%
% First, we define our distance measure

distancemeasure.name = 'ZCR';
distancemeasure.samplerate = 8000;
distancemeasure.mono = true;
distancemeasure.options.nIterations = 0;
distancemeasure.options.nMixtures = 10;
distancemeasure.options.silent = true;

%%
% The 'computefeature' field is called using 'feval'. The wave signal is
% in the variable 'wav'. If an 'options' field was specified, it will be
% transfered to the variable 'options'. The result is returned in the
% variable 'feature'.
distancemeasure.computefeature = 'feature = isp_gmmtrain((wav.^2)'', options);';
%%
% The 'computedistance' field works in the same way. The variables
% 'feature1' and 'feature2' contain two features, and the distance between
% them is returned in the variable 'featureDistance'. Usually,
% 'computefeature' and 'computedistance' will just contain function calls
% and not the entire feature and distance computation code as here.
distancemeasure.computedistance = ['featureDistance=isp_gmmdistance(feature1, ' ...
                    'feature2, options);'];

%%
% Next, we extract the feature from two songs. The isp_extractfeature
% function does all necessary conversion to obtain a 8 kHz mono signal as
% specified. Compressed files are decompressed, and MIDI files are
% synthesized.

file1 = fullfile(isp_toolboxpath, 'Loveshadow - The_Acorns. Seedin Time in The Oak Room - excerpt.mp3');
features{1} = isp_extractfeature(file1, distancemeasure);

file2 = fullfile(isp_toolboxpath, 'shortmidifiles', 'Chill.mid');
features{2} = isp_extractfeature(file2, distancemeasure);

file3 = fullfile(isp_toolboxpath, 'shortmidifiles', 'Latin Five.mid');
features{3} = isp_extractfeature(file3, distancemeasure);

%%
% Finally, we compute the distance matrix and plot it. If the specified
% distance measure has a 'computedistancematrix' field, this is used
% to compute the distance matrix. If not, the method defined in the
% 'computedistance' field is called to compute all entries of the distance
% matrix

dstMat = isp_computedistance(distancemeasure, features, features)

imagesc(dstMat)
colormap(gray)

%%
% If we wanted to store the features in e.g. a database and it is not of
% a simple type such as a vector or string, we can convert it to a vector
% of doubles using the isp_pickle function.

feature1vec = isp_pickle(features{1});

%%
% We can decode it to the original feature:
feature1 = isp_unpickle(feature1vec);

isequalwithequalnans(feature1, features{1})
whos feature1vec feature1
