function [o,options]=isp_mfccvb(s,options)
%ISP_MFCCVB  MFCC implementation based on Mike Brookes' Voicebox
%
% SYNTAX
%   [mfccs,options]=isp_mfccvb(s,options ...)
%
% DESCRIPTION
%   Wrapper function for MELCEPST from Mike Brookes' Voicebox toolbox.
%
% INPUT
%   s:
%     Sound signal.
%   options ...:
%     Field/value pairs or structs with the following fields:
%     mfccprsec:
%       Number of MFCCs per second (default 100).
%     samplerate:
%       Minimum 22050 Hz (default 44100 Hz).
%     nc:
%       Number of mfccs besides the first (default 6).
%     frame:
%       Frame size ( default 20 ms = floor((FS/mfccprsec)*2) ).
%     hopsize:
%       Hop size (default 10 ms = floor(FS/mfccprsec)).
%     nMelbanks:
%       Number of mel banks (default 30).
%     lowercut:
%       Lowest frequency in fraction of sampling rate(default 0).
%     highcut:
%       Higest frequency in fraction of sampling rate (default 11025/FS -> cutoff 55125 hz ).
%     w:
%       Any sensible combination of the following: (default 0Mta)
%       'R'  rectangular window in time domain
%       'N'     Hanning window in time domain
%       'M'     Hamming window in time domain (default)
%
%       't'  triangular shaped filters in mel domain (default)
%       'n'  hanning shaped filters in mel domain
%       'm'  hamming shaped filters in mel domain
%
%       'p'     filters act in the power domain
%       'a'     filters act in the absolute magnitude domain (default)
%
%       '0'  include 0'th order cepstral coefficient
%       'e'  include log energy
%       'd'     include delta coefficients (dc/dt)
%       'D'     include delta-delta coefficients (d^2c/dt^2)
%
%       'z'  highest and lowest filters taper down to zero (default)
%       'y'  lowest filter remains at 1 down to 0 frequency and highest filter remains at 1 up to nyquist freqency
%       If 'ty' or 'ny' is specified, the total power in the fft is preserved.
%
% OUTPUT
%   mfccs:
%     Matrix where mfccs(:,n) is the n'th MFCC vector.
%   options:
%     The sum of input options and default options.
%
% EXAMPLE
%     [mfccs]=isp_mfcc(s)
%     [mfccs,options]=isp_mfcc(s,options)
%
% HISTORY
%   copyrights Intelligent sound 2006 and Mike Brookes
%   licence GPL?
%   author Tue Lehn-Schiøler, ISP,IMM,DTU
%   date 08-02-2006 
%   version 1.0

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


if nargin < 2
    par='';
end
if nargin == 2
    par=options;
end


par=setpars(par);

if par.samplerate<22050;
    warning('Sampling rate below 22050, output cannot be compared to other mfccs')
    par.highcut=0.5;
end

[o]=isp_melcepst(s,par.samplerate,par.w,par.nc,par.nMelbanks,par.frame,par.hopsize,par.lowercut,par.highcut)';

options=par;


function par=setpars(par)
if (~isfield(par,'mfccprsec'))
    par.mfccprsec=100;
end
if (~isfield(par,'samplerate'))
    par.samplerate=44100;
end
FS=par.samplerate;
if (~isfield(par,'nc'))
    par.nc=6;
end
if (~isfield(par,'frame'))
    par.frame=floor((FS/par.mfccprsec)*2);
end
if (~isfield(par,'hopsize'))
    par.hopsize=floor(FS/par.mfccprsec);
end
if (~isfield(par,'nMelbanks'))
    par.nMelbanks=30;
end
if (~isfield(par,'lowercut'))
    par.lowercut=0;
end
if (~isfield(par,'highcut'))
    par.highcut=11025/FS;
end
if (~isfield(par,'w'))
    par.w='0Mta';
end
par.timestamp=datestr(now);


function f=isp_enframe(x,win,inc)
%ENFRAME split signal up into (overlapping) frames: one per row. F=(X,WIN,INC)
%
%       F = ENFRAME(X,LEN) splits the vector X up into
%       frames. Each frame is of length LEN and occupies
%       one row of the output matrix. The last few frames of X
%       will be ignored if its length is not divisible by LEN.
%       It is an error if X is shorter than LEN.
%
%       F = ENFRAME(X,LEN,INC) has frames beginning at increments of INC
%       The centre of frame I is X((I-1)*INC+(LEN+1)/2) for I=1,2,...
%       The number of frames is fix((length(X)-LEN+INC)/INC)
%
%       F = ENFRAME(X,WINDOW) or ENFRAME(X,WINDOW,INC) multiplies
%       each frame by WINDOW(:)

%       Copyright (C) Mike Brookes 1997
%
%      Last modified Tue May 12 13:42:01 1998
%
%   VOICEBOX home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx=length(x);
nwin=length(win);
if (nwin == 1)
   len = win;
else
   len = nwin;
end
if (nargin < 3)
   inc = len;
end
nf = fix((nx-len+inc)/inc);
f=zeros(nf,len);
indf= inc*(0:(nf-1)).';
inds = (1:len);
f(:) = x(indf(:,ones(1,len))+inds(ones(nf,1),:));
if (nwin > 1)
    w = win(:)';
    f = f .* w(ones(nf,1),:);
end


function [x,mn,mx]=isp_melbankm(p,n,fs,fl,fh,w)
%MELBANKM determine matrix for a mel-spaced filterbank [X,MN,MX]=(P,N,FS,FL,FH,W)
%
% Inputs:       p   number of filters in filterbank
%               n   length of fft
%               fs  sample rate in Hz
%               fl  low end of the lowest filter as a fraction of fs (default = 0)
%               fh  high end of highest filter as a fraction of fs (default = 0.5)
%               w   any sensible combination of the following:
%                     't'  triangular shaped filters in mel domain (default)
%                     'n'  hanning shaped filters in mel domain
%                     'm'  hamming shaped filters in mel domain
%
%                     'z'  highest and lowest filters taper down to zero (default)
%                     'y'  lowest filter remains at 1 down to 0 frequency and
%                          highest filter remains at 1 up to nyquist freqency
%
%                      If 'ty' or 'ny' is specified, the total power in the fft is preserved.
%
% Outputs:      x     a sparse matrix containing the filterbank amplitudes
%                     If x is the only output argument then size(x)=[p,1+floor(n/2)]
%                     otherwise size(x)=[p,mx-mn+1]
%               mn   the lowest fft bin with a non-zero coefficient
%               mx   the highest fft bin with a non-zero coefficient
%
% Usage:        f=fft(s);                       f=fft(s);
%               x=melbankm(p,n,fs);             [x,na,nb]=melbankm(p,n,fs);
%               n2=1+floor(n/2);                z=log(x*(f(na:nb)).*conj(f(na:nb)));
%               z=log(x*abs(f(1:n2)).^2);
%               c=dct(z); c(1)=[];
%
% To plot filterbanks e.g.      plot(melbankm(20,256,8000)')
%


%      Copyright (C) Mike Brookes 1997
%
%      Last modified Tue May 12 16:15:28 1998
%
%   VOICEBOX home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=30;
f1=0;
fh=11025/fs;
if fs<22100
    warning('Samling rate below 22100, results may be wrong');
    fh=0.5;
end

if nargin < 6
  w='tz';
  if nargin < 5
    fh=0.5;
    if nargin < 4
      fl=0;
    end
  end
end
f0=700/fs;
fn2=floor(n/2);
lr=log((f0+fh)/(f0+fl))/(p+1);
% convert to fft bin numbers with 0 for DC term
bl=n*((f0+fl)*exp([0 1 p p+1]*lr)-f0);
b2=ceil(bl(2));
b3=floor(bl(3));
if any(w=='y')
  pf=log((f0+(b2:b3)/n)/(f0+fl))/lr;
  fp=floor(pf);
  r=[ones(1,b2) fp fp+1 p*ones(1,fn2-b3)];
  c=[1:b3+1 b2+1:fn2+1];
  v=2*[0.5 ones(1,b2-1) 1-pf+fp pf-fp ones(1,fn2-b3-1) 0.5];
  mn=1;
  mx=fn2+1;
else
  b1=floor(bl(1))+1;
  b4=min(fn2,ceil(bl(4)))-1;
  pf=log((f0+(b1:b4)/n)/(f0+fl))/lr;
  fp=floor(pf);
  pm=pf-fp;
  k2=b2-b1+1;
  k3=b3-b1+1;
  k4=b4-b1+1;
  r=[fp(k2:k4) 1+fp(1:k3)];
  c=[k2:k4 1:k3];
  v=2*[1-pm(k2:k4) pm(1:k3)];
  mn=b1+1;
  mx=b4+1;
end
if any(w=='n')
  v=1-cos(v*pi/2);
elseif any(w=='m')
  v=1-0.92/1.08*cos(v*pi/2);
end
if nargout > 1
  x=sparse(r,c,v);
else
  x=sparse(r,c+mn-1,v,p,1+fn2);
end

function c=isp_melcepst(s,fs,w,nc,p,n,inc,fl,fh)
%MELCEPST Calculate the mel cepstrum of a signal C=(S,FS,W,NC,P,N,INC,FL,FH)
%
%
% Simple use: c=melcepst(s,fs)  % calculate mel cepstrum with 12 coefs, 256 sample frames
%                                 c=melcepst(s,fs,'e0dD') % include log energy, 0th cepstral coef, delta and delta-delta coefs
%
% Inputs:
%     s  speech signal
%     fs  sample rate in Hz (default 11025)
%     nc  number of cepstral coefficients excluding 0'th coefficient (default 12)
%     n   length of frame (default power of 2 <30 ms))
%     p   number of filters in filterbank (default floor(3*log(fs)) )
%     inc frame increment (default n/2)
%     fl  low end of the lowest filter as a fraction of fs (default = 0)
%     fh  high end of highest filter as a fraction of fs (default = 0.5)
%
%               w   any sensible combination of the following:
%
%                               'R'  rectangular window in time domain
%                               'N'     Hanning window in time domain
%                               'M'     Hamming window in time domain (default)
%
%                     't'  triangular shaped filters in mel domain (default)
%                     'n'  hanning shaped filters in mel domain
%                     'm'  hamming shaped filters in mel domain
%
%                               'p'     filters act in the power domain
%                               'a'     filters act in the absolute magnitude domain (default)
%
%                          '0'  include 0'th order cepstral coefficient
%                               'e'  include log energy
%                               'd'     include delta coefficients (dc/dt)
%                               'D'     include delta-delta coefficients (d^2c/dt^2)
%
%                     'z'  highest and lowest filters taper down to zero (default)
%                     'y'  lowest filter remains at 1 down to 0 frequency and
%                                 highest filter remains at 1 up to nyquist freqency
%
%                      If 'ty' or 'ny' is specified, the total power in the fft is preserved.
%
% Outputs:      c     mel cepstrum output: one frame per row. Log energy, if requested, is the
%                 first element of each row followed by the delta and then the delta-delta
%                 coefficients.
%


%      Copyright (C) Mike Brookes 1997
%
%      Last modified Mon May 20 10:35:32 2002
%
%   VOICEBOX is a MATLAB toolbox for speech processing. Home page is at
%   http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 fs=11025; end
if nargin<3 w='M'; end
if nargin<4 nc=12; end
if nargin<5 p=floor(3*log(fs)); end
if nargin<6 n=pow2(floor(log2(0.03*fs))); end
if nargin<9
   fh=0.5;   
   if nargin<8
     fl=0;
     if nargin<7
        inc=floor(n/2);
     end
  end
end

if length(w)==0
   w='M';
end
if any(w=='R')
   z=isp_enframe(s,n,inc);
elseif any (w=='N')
   z=isp_enframe(s,hanning(n),inc);
else
   z=isp_enframe(s,hamming(n),inc);
end
f=isp_rfft(z.');
[m,a,b]=isp_melbankm(p,n,fs,fl,fh,w);
pw=f(a:b,:).*conj(f(a:b,:));
pth=1e-8;
if any(w=='p')
   y=log(max(m*pw,pth));
else
   ath=sqrt(pth);
   y=log(max(m*abs(f(a:b,:)),ath));
end
c=isp_rdct(y).';
nf=size(c,1);
nc=nc+1;
if p>nc
   c(:,nc+1:end)=[];
elseif p<nc
   c=[c zeros(nf,nc-p)];
end
if ~any(w=='0')
   c(:,1)=[];
   nc=nc-1;
end
if any(w=='e')
   c=[log(sum(pw)).' c];
   nc=nc+1;
end

% calculate derivative

if any(w=='D')
  vf=(4:-1:-4)/60;
  af=(1:-1:-1)/2;
  ww=ones(5,1);
  cx=[c(ww,:); c; c(nf*ww,:)];
  vx=reshape(filter(vf,1,cx(:)),nf+10,nc);
  vx(1:8,:)=[];
  ax=reshape(filter(af,1,vx(:)),nf+2,nc);
  ax(1:2,:)=[];
  vx([1 nf+2],:)=[];
  if any(w=='d')
     c=[c vx ax];
  else
     c=[c ax];
  end
elseif any(w=='d')
  vf=(4:-1:-4)/60;
  ww=ones(4,1);
  cx=[c(ww,:); c; c(nf*ww,:)];
  vx=reshape(filter(vf,1,cx(:)),nf+8,nc);
  vx(1:8,:)=[];
  c=[c vx];
end
 
if nargout<1
   [nf,nc]=size(c);
   t=((0:nf-1)*inc+(n-1)/2)/fs;
   ci=(1:nc)-any(w=='0')-any(w=='e');
   imh = imagesc(t,ci,c.');
   axis('xy');
   xlabel('Time (s)');
   ylabel('Mel-cepstrum coefficient');
        map = (0:63)'/63;
        colormap([map map map]);
        colorbar;
end


function y=isp_rdct(x,n)
%RDCT     Discrete cosine transform of real data Y=(X,N)
% Data is truncated/padded to length N.
%
% This routine is equivalent to multiplying by the matrix
%
%   rdct(eye(n)) = diag([sqrt(2) 2*ones(1,n-1)]) * cos((0:n-1)'*(0.5:n)*pi/n)
%
% The rows and columns of the matrix are orthogonal but not unit modulus.
% Various versions of the DCT are obtained by pre-multiplying the above
% matrix by diag([b/a ones(1,n-1)/a]) and post-multiplying the
% inverse transform matrix by its inverse. A common choice is a=n and/or b=sqrt(2).
% Choose a=sqrt(2n) and b=1 to make the matrix orthogonal.
% If b~=1 then the columns are no longer orthogonal.
%
% see IRDCT for the inverse transform




%      Copyright (C) Mike Brookes 1998
%
%      Last modified Tue Apr 13 15:56:48 1999
%
%   VOICEBOX is a MATLAB toolbox for speech processing. Home page is at
%   http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fl=size(x,1)==1;
if fl x=x(:); end
[m,k]=size(x);
if nargin<2 n=m;
elseif n>m x=[x; zeros(n-m,k)];
elseif n<m x(n+1:m,:)=[];
end

x=[x(1:2:n,:); x(2*fix(n/2):-2:2,:)];
z=[sqrt(2) 2*exp((-0.5i*pi/n)*(1:n-1))].';
y=real(fft(x).*z(:,ones(1,k)));

if fl y=y.'; end

function y=isp_rfft(x,n,d)
%RFFT     FFT of real data Y=(X,N)
% Data is truncated/padded to length N if specified.
%   N even:     (N+2)/2 points are returned with
%                       the first and last being real
%   N odd:      (N+1)/2 points are returned with the
%                       first being real



%      Copyright (C) Mike Brookes 1998
%
%      Last modified Fri Mar  7 15:43:06 2003
%
%   VOICEBOX is a MATLAB toolbox for speech processing. Home page is at
%   http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=size(x);
if prod(s)==1
    y=x
else
    if nargin <3
        d=find(s>1);
        d=d(1);
        if nargin<2
            n=[];
        end
    end
    y=fft(x,n,d);
    y=reshape(y,prod(s(1:d-1)),s(d),prod(s(d+1:end))); 
    s(d)=1+fix(s(d)/2);
    y(:,s(d)+1:end,:)=[];
    y=reshape(y,s);
end

