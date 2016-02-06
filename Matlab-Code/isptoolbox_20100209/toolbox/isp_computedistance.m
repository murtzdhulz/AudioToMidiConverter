%ISP_COMPUTEDISTANCE  Compute distance matrix.
%
% SYNTAX
%   [distancematrix, comptime] = isp_computedistance(distance, features1, features2)
%
% DESCRIPTION
%   Computes the distance matrix from features previously extracted.
%
% INPUT
%   distance:
%     Struct defining a distance measure as returned by
%     e.g. isp_mfccgmmkl. This function uses the 'options' and either the
%     'computedistancematrix' or 'computedistance' fields.
%   features1, features2:
%     Cell array of features.
%
% OUTPUT
%   distancematrix:
%     Matrix where entry (m, n) is the distance between features1{m} and
%     features2{n}
%   comptime:
%     Computation time.
%
% SEE ALSO
%   isp_evaluate, isp_extractfeature.
%
% HISTORY
%   Created by Jesper Højvang Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.

function [distancematrix, distancematrixcomputationtime]=isp_computedistance(distance, features1, features2)
    
    if ~iscell(features1)
        features1={features1};
    end
    if ~iscell(features2)
        features2={features2};
    end

    if isfield(distance, 'options')
        options = distance.options;
    else
        options = struct;
    end

    if isfield(distance, 'computedistancematrix')
        cputimestart = cputime;
        eval(distance.computedistancematrix);
        distancematrixcomputationtime = cputime - cputimestart;
    else
        distancecmd = distance.computedistance;
        distancematrix = zeros(length(features1), length(features2));
        cputimestart = cputime;
        for iFeature1=1:length(features1)
            for iFeature2=1:length(features2)
                feature1 = features1{iFeature1};
                feature2 = features2{iFeature2};
                eval(distancecmd);
                distancematrix(iFeature1, iFeature2) = featureDistance;
            end
        end
        distancematrixcomputationtime = cputime - cputimestart;
    end

end
