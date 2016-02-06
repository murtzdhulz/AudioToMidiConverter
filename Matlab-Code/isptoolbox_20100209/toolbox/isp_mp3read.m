%ISP_MP3READ  Read wave data and song information from an MP3 file. 
%
% SYNTAX
%   [y,info,options] = isp_mp3read(mp3file,options)
%
% DESCRIPTION
%   This script reads an mp3 file, decode it, and return the
%   wav signal and optionally returns song information from the ID3 tags as
%   well.
%
% INPUT
%   mp3file:
%     Input file name.
%   options ...:
%     Structs or field/value pairs specifying any of the following.
%     downsamp:
%       A number indicating if the content of the mp3 file
%       should be downsampled or not. 1 indicates no downsampling, 2
%       indicates half, and 4 indicates quarter sampling
%       frequency. (Default: 1)
%     forcemono:
%       A number indicating if the content of the mp3 file
%       should be forced to mono or kept as original. 0 indicates no
%       forced mono, while 1 indicates forced mono. (Default: 0)
%     nbits:
%       A number indicating number of bits per sample. (Default: 16)
%     toolboxpath:
%       A textstring indicating the path to the folder holding the 
%       mp3 read tools. (Default: autodetect)
%
% OUTPUT
%   y:
%     Sampled audio.
%   info:
%     Struct with information from ID3 tags.
%   options:
%     Similar to the 'options' input argument, but with added defaults
%     for unspecified values.
%
% EXAMPLE
%     Y = isp_mp3read(mp3file)
%     [Y,info] = isp_mp3read(mp3file)
%     [Y,info,options] = isp_mp3read(mp3file,options)
%
% HISTORY
%   2003-07-20:  dpwe@ee.columbia.edu  This version used mpg123.
%   2004-08-31:  Fixed to read whole files correctly
%   2004-09-08:  Uses mp3info to get info about mp3 files too
%   2004-09-18:  Reports all mp3info fields in OPTS.fmt; handles MPG2LSF sizes
%                + added MONO, DOWNSAMP flags, changed default behavior.
%   2006-02-13:  Strongly modified by Kaare Brandt Petersen. Cleaned and
%                modified to fit ISP Toolbox format and made more robust
%                to errors and deviations. 
%   2006-03-15:  Minor bugfix by tls, default monoing value actually set to 0
%   2006-03-15:  Modified by tls to return year artist etc in the info struct
%   2006-03-30:  Minor modification by jhj: forcemono help text changed, and
%                linux support added.
%   2006-05-08:  Major modification, now uses java code by Elvis Ramonovic to
%                extract mpeginfo.
%

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function [Y,info,options] = isp_mp3read(mp3file,varargin)

if nargin < 2
    options=struct;
end

options=isp_interpretarguments(struct(...
    'toolboxpath', isp_toolboxpath, ...
    'forcemono', 0, ... % 0 = no forced mono (keep original) and 1 = yes, forced mono
    'downsamp', 1, ... % 1 = no downsampling, 2 = half, 4=qtr
    'nbits', 16, ...
    'timestamp', datestr(now)), varargin{:});

% Settings defaults values of variables
toolboxpath = options.toolboxpath;
forcemono   = options.forcemono;
downsamp    = options.downsamp;
nbits       = options.nbits;
info.error = '';


%%%%%% Location of the binaries
switch(computer)
    case 'GLNX86'
        madplay = ['"' toolboxpath filesep 'isp_madplay.glnx86"'];
        deletecommand = 'rm ';
    case 'PCWIN'
        madplay = ['"' toolboxpath filesep 'isp_madplay"'];

        deletecommand = 'del ';
    case 'SOL2'
        madplay = ['"' toolboxpath filesep 'isp_madplay.sol2"'];
        deletecommand = 'rm ';
    otherwise
        error(['Sorry, the platform ' computer ' is not supported.'])
end

%Add java to class path
javaaddpath([toolboxpath filesep 'MatlabSongInfo.jar']);

try 
    k = MPEGinfo(mp3file);
catch
    [pdummy,ndummy,edummy,vdummy] = fileparts(mp3file);
    mystr = lasterr; 
    warning([mystr '\n'])
    Y = 0;
    info.error = mystr;
    return
end

if (~k.getStatus)
    vv=char(k.getErrorMessage);
   warning(vv)
end

SR = k.getFrequency;
nframes = k.getFrames;

%k.getMode;
%nchans = 2 - strcmp(val1{9}, 'mono');
nchans=2;

layer = char(k.getLayerVersion);
% bitrate = nums{3}*1000;
mpgv = k.getMPEGVersion;

% Figure samples per frame, after
% http://board.mp3-tech.org/view.php3?bn=agora_mp3techorg&key=1019510889
if layer == 1
    smpspfrm = 384;
elseif SR < 32000 & layer ==3
    smpspfrm = 576;
    if mpgv == 1
        error('SR < 32000 but mpeg version = 1');
    end
else
    smpspfrm = 1152;
end


% OPTS.fmt.mpgBitrate = bitrate;
OPTS.fmt.mpgVersion = mpgv;
% fields from wavread's OPTS
OPTS.fmt.nSamplesPerSec = SR;
OPTS.fmt.nChannels = nchans;
OPTS.fmt.nBitsPerSample = nbits;
OPTS.fmt.mpgNFrames = nframes;
OPTS.fmt.mpgCopyright = k.getCopyrighted;
%OPTS.fmt.mpgEmphasis = val1{5};
OPTS.fmt.mpgCRC =k.getCrcs;
OPTS.fmt.mpgLayer = k.getLayerVersion;
OPTS.fmt.mpgOriginal = k.getOriginal;
OPTS.fmt.mpgChanmode = k.getMode;
%OPTS.fmt.mpgPad = val1{10};
OPTS.fmt.mpgSampsPerFrame = smpspfrm;

if SR == 16000 & downsamp == 4
    error('mpg123 will not downsample 16 kHz files by 4 (only 2)');
end

downsampstr = [' -',num2str(downsamp)];
FS = SR/downsamp;

if forcemono == 1
    nchans = 1;
    chansstr = ' -m';
else
    chansstr = '-S';
    nchans = 2;
end



N(1)=1;
N(2)=floor(smpspfrm*nframes/downsamp);



% Temporary file to use
tmpdir = '';
tmpfile = [tmpdir num2str(round(1000*rand(1))) '.wav'];

skipx = 0;
skipblks = 0;
skipstr = '';
sttfrm = N(1)-1;
if sttfrm > 0
    skipblks = floor(sttfrm*downsamp/smpspfrm);
    skipx = sttfrm - (skipblks*smpspfrm/downsamp);
    skipstr = [' -k ', num2str(skipblks)];
end



ssec = rem( fix(N(1)/FS) ,  60) ;   % start time seconds
smin = fix( (N(1)/FS) / 60 );    % start time minutes
durtime = fix(N(2)/FS) - fix(N(1)/FS) ; % total duration in seconds
dursec = rem( durtime, 60 ); % duration
durmin = fix( durtime / 60 );
if smin < 10,
    mytime{1} = ['0' num2str(smin)];
else
    mytime{1} = num2str(smin);
end
if ssec < 10,
    mytime{2} = ['0' num2str(ssec)];
else
    mytime{2} = num2str(ssec);
end
if durmin < 10,
    mytime{3} = ['0' num2str(durmin)];
else
    mytime{3} = num2str(durmin);
end
if dursec < 10,
    mytime{4} = ['0' num2str(dursec)];
else
    mytime{4} = num2str(dursec);
end




% Run the decode
% cmd=[mpg123, downsampstr, chansstr, skipstr, lenstr, ' -q -w ', tmpfile, ' ', '"',mp3file,'"']
% NOTE: Madplay seemed to do the job - contrary to mpg123.. at least on a win32 system
%cmd=[madplay ' -S ' ' --start=00:' mytime{1} ':' mytime{2} '  --time=00:' mytime{3} ':' mytime{4} ' --output="wave:' tmpfile '" ' '"' mp3file '"'];
%disp(cmd);
%chansstr = '-S';
cmd=[madplay ' -2 ' chansstr ' --start=00:' mytime{1} ':' mytime{2} '  --time=00:' mytime{3} ':' mytime{4} ' --output="wave:' tmpfile '" ' '"' mp3file '"'];





try 
    mysystem(cmd);
catch
    [pdummy,ndummy,edummy,vdummy] = fileparts(mp3file);
    mystr = lasterr; %sprintf('ERROR: isp_mp3read has terminated processing "%s%s" because of a serious problem with "madplay.exe".',ndummy,edummy);
    warning([mystr '\n'])
    Y = 0;
    info.error = mystr;
    return
end

% Load the data
[Y,SR] = wavread(tmpfile);

% Delete tmp file
mysystem([deletecommand tmpfile]);


%1 %Q 2 %u 3 %v 4 %C 5 %e 6 %E 7 %L 8 %O 9 %o 10 %p 11 %G 12 %n 13 %t 14 %S
%15 %l"

%info.length = k.getLength;
info.samplingfrequency = FS;
info.nbits = nbits;
info.other = OPTS.fmt;

if (~k.tags.getId3v1);
    info.id3v1.status=0;
else
    info.id3v1.status=1;
    % info.genrenum = num2str(val1{11});
    info.id3v1.tracknumber = char(k.tags.getTrackID3v1);
    info.id3v1.title = char(k.tags.getTitleID3v1);

    info.id3v1.album =  char(k.tags.getAlbumID3v1);
    info.id3v1.comment = char(k.tags.getCommentID3v1);
    info.id3v1.year = char(k.tags.getYearID3v1);
    info.id3v1.author = char(k.tags.getArtistID3v1);
    info.id3v1.genre = char(k.tags.getGenreID3v1);
end

%if (k.tags.getID3v2);
if(0)
    info.id3v2.status=0;
else
    info.id3v2.status=1;
    % info.genrenum = num2str(val1{11});
    info.id3v2.tracknumber = char(k.tags.getTrackID3v2);
    info.id3v2.title = char(k.tags.getTitleID3v2);

    info.id3v2.album =  char(k.tags.getAlbumID3v2);
    info.id3v2.comment = char(k.tags.getCommentID3v2);
    info.id3v2.year = char(k.tags.getYearID3v2);
    info.id3v2.author = char(k.tags.getArtistID3v2);
    info.id3v2.genre = char(k.tags.getGenreID3v2);
    info.id3v2.composer=char(k.tags.getComposerID3v2);
    info.id3v2.copyright=char(k.tags.getCopyRigthID3v2);
    info.id3v2.encodedby=char(k.tags.getEncodedByID3v2);
    info.id3v2.originalartist=char(k.tags.getOrrgArtistID3v2);
    info.id3v2.url=char(k.tags.getUrlID3v2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = mysystem(cmd)
% Run system command; report error; strip all but last line
[s,w] = system(cmd);
if s ~= 0
    error(['unable to execute ',cmd]);
end
% Keep just final line
w = w((1+max([0,findstr(w,10)])):end);
% Debug
%disp([cmd,' -> ','*',w,'*']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = tokenize(s)
% Break space-separated string into cell array of strings
% 2004-09-18 dpwe@ee.columbia.edu
a = [];
p = 1;
n = 1;
l = length(s);
nss = findstr([s(p:end),' '],' ');
for ns = nss
    % Skip initial spaces
    if ns == p
        p = p+1;
    else
        if p <= l
            a{n} = s(p:(ns-1));
            n = n+1;
            p = ns+1;
        end
    end
end

