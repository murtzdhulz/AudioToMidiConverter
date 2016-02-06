%ISP_ARWRAPPER  Wrapper function for isp_ar.
%
% SYNTAX
%   [AR,ceps,options]=isp_arwrapper(s, options ...);
%
% DESCRIPTION
%   Extracts MFCCs and then AR cooeficients.
%
% INPUT
%   s:
%     Wave signal.
%   options ...:
%     Field/value pairs or structs with the following field names:
%     ar:
%       Struct specifying options for isp_ar.
%     mfcc:
%       Struct specifying options for isp_mfccvb.
%
% OUTPUT
%   AR:
%     Matrix vith AR coeeficients.
%   ceps:
%     Matrix with MFCCs.
%   options:
%     Input options plus default options.
%
% EXAMPLE
%     mfccopts.samplerate = 22050;
%     [Ar,mfcc,options]=isp_arwrapper(rand(1,100000), 'mfcc', mfccopts);
%
% SEE ALSO
%   isp_ar, isp_mfccvb.
%
% HISTORY
%   Copyrights Intelligent sound 2006
%   Author Tue Lehn-Schiøler, ISP,IMM,DTU
%   date 13-02-2006
%
%   Modified 15/03/06 by TLS, Automatically convert stereo to mono
%   Modified by Jesper Højvang Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [AR,ceps,options]=isp_arwrapper(s,varargin);

options = isp_interpretarguments(struct('ar', struct, 'mfcc', struct), varargin{:});

if size(s,2)==2
    s=mean(s,2);
end

% extract mfccs with standard settings
[ceps,options.mfcc]=isp_mfccvb(s,options.mfcc);
%[ceps]=mfcc_ar(Y+10*eps,FS,mfccprsec,20);

mfccprsec=options.mfcc.mfccprsec;

nn=length(ceps);
n=(nn-mod(nn,mfccprsec))/mfccprsec;
AR=[];

if(~isfield(options,'ar'))
    options.ar='';
end
if(~isfield(options.ar,'p'))
    options.ar.p=3;
end
%remove first mfcc due to power issues
ceps2=ceps(:,2:end);
for i=1:n
    idxx=(i-1)*mfccprsec+1:i*mfccprsec;
    tmpmfcc=ceps2(idxx,:)';
    [A,w,C,res,options] = isp_ar(tmpmfcc,options.ar);
    ARtemp = isp_ar2vector(A,w,C);
%     c=nonzeros(triu(C));
%     ARtemp = [A(:);w(:);c];
    AR=[AR ARtemp];
end


function Arvec=isp_ar2vector(A,w,C);
% @function 
% Arvec=isp_ar2vector(A,w,C);
%
%@description
% returns the AR coefficients as a vector
%
%%@input 
% X : (D x N) multidimensional times-series input
%
%@output
% Arvec: vector with ar cooeficients
%
%@author Anders Meng, IMM, DTU, DK.
%@date 27-05-05 
%@version 1.0
%
%@changelog
% modified to fit toolbox format by Tue Lehn-Schiøler 13/02-06
% modified to also work with silent input by Jesper H. Jensen
%
%@examples
% Arvec=isp_ar2vector(A,w,C);

c=nonzeros(triu(C));
Arvec = [A(:);w(:);c];