%isp_extractfeature  Extract features from a song.
%
% SYNTAX
%   [feature, song, modwav, options] = isp_extractfeature(song, distancemeasure, options ...)
%
% DESCRIPTION
%   Extracts the specified feature from a song
%
% INPUT
%   song:
%     Either the filename of a song or a struct as returned by
%     isp_makesonglist.
%   distancemeasure:
%     A struct defining a distance measure as returned by
%     e.g. isp_mfccgmmkl. The fields mono, samplerate, computefeature and
%     options are used.
%   options ...:
%     Structs or field/value pairs with any of the following field names:
%     addsilence, normalize, transpose, duration, or percussion:
%       If the input song is a MIDI file, isp_midimodify is called
%       before extracting the feature.
%     removesilence, bitrate, snr, bandwidth, maxlength:
%       Modify the audio data by calling isp_wavmodify with the
%       specified options before extracting the feature.
%     soundfont:
%       The sound font used to synthesize MIDI files.
%     maxSamples:
%       This option is passed on to isp_audioread. Unlike the maxlength
%       option, this option takes effect before the wave data is read
%       into memory.
%
% OUTPUT
%   feature:
%     The extracted feature.
%   song:
%     Struct specifying the modifications applied before feature extraction.
%   modwav:
%     Raw audio data from which the feature was extracted.
%   options:
%     Input options supplemented with default values.
%
% SEE ALSO
%   isp_audioread, isp_computedistance, isp_midimodify, isp_wavmodify.
%
% HISTORY
%   Created by Jesper Højvang Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.

function [feature, song, modwav, options] = isp_extractfeature(song, distancemeasure, varargin)

    if ischar(song)
        cellout = false;
        temp = song;
        song = struct;
        song.filename = temp;
        song.modification = struct;
        song.soundfont = '';
    else
        cellout=true;
    end

    if ~isfield(distancemeasure, 'options')
        distancemeasure.options = struct;
    end

    modopts = struct('addsilence', 0, ...  % Midi options
                     'normalize', [], ...
                     'transpose', 0, ...
                     'duration', 1, ...
                     'percussion', false, ...
                     'removesilence', true, ... % Wav options
                     'bitrate', inf, ...
                     'snr', inf, ...
                     'bandwidth', inf, ...
                     'maxlength', inf, ...
                     'soundfont', '', ...
                     'maxSamples', inf);

    options = isp_interpretarguments(modopts, varargin{:});
    
    % Order of precedence: function arguments, song.modification, defaults from modopts.

    prevFilename='';
    prevMidi = [];
    prevSf = '';


    if ~isfield(song(1), 'modification')
        for iSong=1:numel(song)
            song(iSong).modification = struct;
        end
    end

    feature = cell(size(song));
    featureextractiontime = zeros(size(song));
    for iSong=1:numel(song)
        if isfield(song(iSong).modification, 'soundfont')
            warning('Songs shouldn''t have a modification.soundfont field.')
        end

        song(iSong).modification=isp_interpretarguments(modopts, song(iSong).modification, varargin{:});
        song(iSong).modification.mono = distancemeasure.mono;
        if ~isempty(song(iSong).modification.soundfont)
            song(iSong).soundfont = song(iSong).modification.soundfont;
        end
        song(iSong).modification = rmfield(song(iSong).modification, ...
                                           {'soundfont', 'maxSamples'});

        [dummy, dummy, ext] = fileparts(song(iSong).filename);
        switch ext
          case {'.mid'}
            % MIDI
            if strcmp(song(iSong).filename, prevFilename)
                fprintf(1, 'Reusing midi file.\n');
            else
                prevFilename=song(iSong).filename;
                midi=isp_midiread(song(iSong).filename);
            end
            [modmidi, song(iSong).modification]=isp_midimodify(midi, song(iSong).modification);
            if isfield(song(iSong), 'instrument')
                fprintf(1, 'Modifying instrumentation.\n')
                modmidi.instruments(modmidi.instruments(:, 2)~=10, 3) = song(iSong).instrument;
            end
            if isequal(modmidi, prevMidi) && ...
                    isequal(song(iSong).soundfont, prevSf)
                fprintf(1, 'Reusing wav data.\n');
            else
                prevMidi = modmidi;
                prevSf = song(iSong).soundfont;

                if isfield(distancemeasure, 'separatechannels') ...
                        && distancemeasure.separatechannels
                   % Play each channel separately
                   channels = unique(modmidi.notes(:, 3));
                   wav={};
                   for iCh = 1:length(channels)
                       ch = channels(iCh);
                       chNotes = modmidi.notes(:, 3)==ch;
                       tmpmid = modmidi;
                       tmpmid.notes(~chNotes,:) = [];
                       [wav{iCh}, fs]=isp_midisynth(tmpmid, ...
                           song(iSong).soundfont, ...
                           'samplerate', distancemeasure.samplerate, ...
                           'nChannels', 2);
                   end
                
                else                
                    [wav, fs]=isp_midisynth(modmidi, song(iSong).soundfont, ...
                              'samplerate', distancemeasure.samplerate, ...
                              'nChannels', 2);
                end
            end

          otherwise
            % Assume it's sampled audio
            [wav, fs]=isp_audioread(song(iSong).filename, ...
                                    'samplerate', distancemeasure.samplerate, ...
                                    'nChannels', 2, 'maxSamples', options.maxSamples);
        end

        if fs ~= distancemeasure.samplerate
            error(['Sorry, but sampling rates don''t match. This ' ...
                   'is not supposed to happen ...'])
        end

        if iscell(wav)
            modwav={};
            for n=1:numel(wav)
                [modwav{n}, song(iSong).modification]=isp_wavmodify(wav{n}, ...
                    song(iSong).modification, 'samplerate', fs);
            end
        else
            [modwav, song(iSong).modification]=isp_wavmodify(wav, ...
                song(iSong).modification, 'samplerate', fs);
        end

        [feature{iSong}, featureextractiontime(iSong)] = ...
            internalfunctionextractfeatures(distancemeasure, modwav);

        fprintf(1, '\n');
    end
    
    if ~cellout
        feature=feature{1};
    end
    
end

function [feature, featureextractiontime] = internalfunctionextractfeatures(distancemeasure, wav)
% This is a separate function to minimize the risk of overwriting variables
    featurecmd = distancemeasure.computefeature;
    options = distancemeasure.options;
    var = []; % Free variable for use in featurecmd. MATLAB
              % unfortunately no longer allows dynamic creation of
              % new variables in functions.
    cputimestart = cputime;
    eval(featurecmd);
    featureextractiontime = cputime - cputimestart;
end
