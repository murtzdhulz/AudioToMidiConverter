%% ISP_AUDIOREADDEMO  Demonstration of how to use isp_audioread
% Specify full file name of the demo file, which is an excerpt from
% Loveshadow: The Acorns. Seedin Time in The Oak Room.
% Get the entire file from http://ccmixter.org/media/files/Loveshadow/12192

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.

mp3file=fullfile(isp_toolboxpath, 'Loveshadow - The_Acorns. Seedin Time in The Oak Room - excerpt.mp3')

%%
% Read it and play a short excerpt
[wav, fs]=isp_audioread(mp3file);
% wavplay only works on the Windows platform
%wavplay(wav, fs) 

%%
% Do it again, but this time convert to 16kHz mono
[wav, fs]=isp_audioread(mp3file, 'nChannels', 1, 'samplerate', 16000);
% wavplay only works on the Windows platform
%wavplay(wav, fs)

%%
% Use isp_mp3read
[wav, fmt]=isp_audioread(mp3file);
fmt
