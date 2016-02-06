%ISP_IFPTRACK  Pitch track based on instantaneous frequency.
%
% SYNTAX
%   [p,m,S] = ifptrack(d,w,sr,fminl,fminu,fmaxl,fmaxu)
% 
% DESCRIPTION
%   Track pitch based on instantaneous frequency. It looks for adjacent
%   bins with same inst freq.
%
% INPUT
%   d:
%     Input waveform.
%   w:
%     STFT DFT length (window is half, hop is 1/4).
%   sr:
%     Sample rate.
%   fminl, fmaxu, fmaxl, fmaxu:
%     Define ramps at edge of sensitivity.
%
% OUTPUT
%   p:
%     Tracked pitch frequencies.
%   m:
%     Tracked pitch magnitudes.
%   S:
%     The underlying complex STFT.
%
% SEE ALSO
%   isp_ifchromagram, isp_ifgram.
%
% HISTORY
%   2006-05-03:  dpwe@ee.columbia.edu
%   2007-10-11:  Optimized by Jesper Højvang Jensen (jhj@es.aau.dk)

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.



% Dan's original copyright:
%
%   Copyright (c) 2006 Columbia University.
% 
%   This file is part of LabROSA-coversongID
% 
%   LabROSA-coversongID is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License version 2 as
%   published by the Free Software Foundation.
% 
%   LabROSA-coversongID is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with LabROSA-coversongID; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
%   02110-1301 USA
% 
%   See the file "COPYING" for the text of the license.

function [p,m,S] = isp_ifptrack(d,w,sr,fminl,fminu,fmaxl,fmaxu)

% downweight fundamentals below here
if nargin < 4; fminl = 150; end
if nargin < 5; fminu = 300; end
% highest frequency we look to
if nargin < 6; fmaxl = 2000; end
if nargin < 7; fmaxu = 4000; end


% Only look at bins up to 2 kHz
maxbin = round(fmaxu * (w/sr) );
minbin = round(fminl * (w/sr) );

% Calculate the inst freq gram
[I,S] = isp_ifgram(d,w,w/2,w/4,sr, maxbin);


% Compile mex file if necessary
mexname = [mfilename '_helper'];
if ~(exist(mexname)==3)
    if exist('OCTAVE_VERSION')
        npname = fullfile(isp_toolboxpath, mexname);
        mkoctfile([npname '.c'], '--mex', '--output', [npname '.' mexext])
    else
        cfile=which([mexname '.c']);
        try
            mex(cfile, '-outdir', fileparts(cfile))
        catch
            warning('Unable to compile mex file')
        end
    end
end


if exist('isp_ifptrack_helper')==3
    [p,m] = isp_ifptrack_helper(I, S, .75*sr/w, fminl,fminu,fmaxl,fmaxu);
else
    fprintf(1, 'Couldn''t find mex-file. Using matlab version for %s.\n', mfilename);
    % Find plateaus in ifgram - stretches where delta IF is < thr
    ddif = I([2:maxbin maxbin], :) - I([1 1:(maxbin-1)],:);

    % expected increment per bin = sr/w, threshold at 3/4 that
    dgood = abs(ddif) < .75*sr/w;

    % delete any single bins (both above and below are zero);
    dgood = dgood .* (dgood([2:maxbin maxbin],:) > 0 | dgood([1 1:(maxbin-1)],:) > 0);

    % reconstruct just pitchy cells?
    %r = istft(p.*S,w,w/2,w/4);

    p = 0*dgood;
    m = 0*dgood;

    % For each frame, extract all harmonic freqs & magnitudes
    lds = size(dgood,1);
    for t = 1:size(I,2)
        ds = dgood(:,t)';
        % find nonzero regions in this vector
        st = find(([0,ds(1:(lds-1))]==0) & (ds > 0));
        en = find((ds > 0) & ([ds(2:lds),0] == 0));
        npks = length(st);
        frqs = zeros(1,npks);
        mags = zeros(1,npks);
        for i = 1:length(st)
            bump = abs(S(st(i):en(i),t));
            mags(i) = sum(bump);
            frqs(i) = (bump'*I(st(i):en(i),t))/(mags(i)+(mags(i)==0));
            if frqs(i) > fmaxu
                mags(i) = 0;
                frqs(i) = 0;
            elseif frqs(i) > fmaxl
                mags(i) = mags(i) * max(0, (fmaxu - frqs(i))/(fmaxu-fmaxl));
            end
            % downweight magnitudes below? 200 Hz
            if frqs(i) < fminl
                mags(i) = 0;
                frqs(i) = 0;
            elseif frqs(i) < fminu
                % 1 octave fade-out
                mags(i) = mags(i) * (frqs(i) - fminl)/(fminu-fminl);
            end
            if frqs(i) < 0 
                mags(i) = 0;
                frqs(i) = 0;
            end
            
        end

        % then just keep the largest at each frame (for now)
        %  [v,ix] = max(mags);
        %  p(t) = frqs(ix);
        %  m(t) = mags(ix);
        % No, keep them all
        %bin = st;
        bin = round((st+en)/2);
        p(bin,t) = frqs;
        m(bin,t) = mags;
    end

    %% Pull out the max in each column
    %[mm,ix] = max(m);
    %% idiom to retrieve different element from each column
    %[nr,nc] = size(p);
    %pp = p((nr*[0:(nc-1)])+ix);
    %mm = m((nr*[0:(nc-1)])+ix);
    % r = synthtrax(pp,mm,sr,w/4);

    %p = pp;
    %m = mm;
end