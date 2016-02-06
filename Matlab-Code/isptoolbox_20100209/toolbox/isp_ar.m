%ISP_AR  Compute multivariate autoregressive coefficients.
%
% SYNTAX
%   [A,w,C,res,options]=isp_ar(X,options ...)
%
% DESCRIPTION
%   For the model
%     y(n) = A_1 y(n-1) + A_2 y(n-2) + ... + A_p y(n-p) + w + u(n)
%
%   In the model y(n) is the time-series and is a vector of dimension (D x
%   1). The matrices (A_1, A_2, .... , A_p) is determined by the estimation
%   procedure. v is known as intercept terms and u(n) is the white noise
%   process of dimension (D x 1) with zero mean and bounded fourth  
%   order moment. 
%   The model order is assumed to be known for the given problem.
%
% INPUT 
%   X:
%     (D x N) multidimensional times-series input
%   options ...:
%     Parameter/value pairs or struct(s) with field names as follows:
%     p:
%       model order (default 1)
%     type:
%       Either 'var' for vector autoregressive model (((DEFAULT))), or
%       'dar' for diagonal autoregressive model (independence among
%       dimensions).
%
% OUTPUT
%   A:
%     (D x Dp) A matrices stacked.
%   w:
%     (D x 1) intercept term of time-series
%   C:
%     (D x D) Covariance estimate of u(n), hence <u(n)u(n)^T>.
%   res:
%     residuals of predictor.
%   options:
%     Struct with the options used, including default values for
%     unspecified options.
%
% EXAMPLE
%   [A,w,C,res]=isp_ar(X)
%
% SEE ALSO
%   isp_arwrapper
%
% HISTORY
%   27-05-05: Written by Anders Meng, IMM, DTU, DK.
%   13/02-06: modified to fit toolbox format by Tue Lehn-Schiøler.
%   Further modified by Jesper H. Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.

function [A,w,C,res,options]=isp_ar(X,varargin)


options = isp_interpretarguments(struct(...
    'p', 2, ...
    'type', 'var', ...
    'timestamp', datestr(now)), varargin{:});


p=options.p;
if(strcmp(options.type,'dar'))
    type = 2;
elseif(strcmp(options.type,'var'))
    type = 1;
else
    warning('Unrecognized type using VAR')
end

%--
%try
[D,N]=size(X);
if(type==1)
    %-- Vector Autoregressive Model - Estimation using Least Squares (or
    %   Generalized Least Squares)

    %-- Init
    Sk = eye(D); 

    %-- Tradiontal LS estimate (Same as maximum likelihood)
    Z = StackMatrixLocal(X,(1:p)');
    [DZ,NZ]=size(Z);
    Z = [ones(1,NZ);Z];
    X = X(:,p+1:end);
    B = X*Z'/(Z*Z'+sqrt(eps)*eye(DZ+1)); %Some regularization, just in case).
    res = (X-B*Z);
    Sk = 1/(N-p*D-1)*res*res'; %Consistent estimator (unbiased)

    %-- Output
    C = Sk; %Unbiased Covariance matrix estimate
    A = B(:,2:end); %A-matrices (Stacked horizontally)
    w = B(:,1); %Intercept terms

elseif(type==2)
    %-- Diagonal Autoregressive Model - Estimation using Least Squares
    %   (Generalized Least Squares).

    
    for dimn = 1:D     
        Z = StackMatrixLocal(X(dimn,:),(1:p)');
        [DZ,NZ]=size(Z);
        Z = [ones(1,NZ);Z];
        b(dimn,:)=((Z*Z'+sqrt(eps)*eye(DZ+1))\(Z*X(dimn,p+1:end)'))';
        res(dimn,:) = X(dimn,p+1:end) - b(dimn,:)*Z;
        C(dimn)=1/(N-p-1) * res(dimn,:)*res(dimn,:)';
    end
    
    A = b(:,2:end);
    w = b(:,1);
else
    error('Doooo')
end



% catch
% disp('Error in wlsar.m - catch statement')    
% end

function [Z]=StackMatrixLocal(X,I)
% usage : Z = StackMatrix(X,I)
%
%	Where I = [1,2,3] of length P 
%	      X = D x N matrix of data
%	      Z = Dp x N-p	
%
% Anders Meng, Feb.06.
%
[D,N]=size(X);
P = length(I);
Z = zeros(P*D,N-P);
for i=1:P
    Z((i-1)*D+1:i*D,:)=X(:,P-i+1:end-i);
end