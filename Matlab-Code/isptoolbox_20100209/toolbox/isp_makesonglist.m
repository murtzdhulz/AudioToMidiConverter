function [songArray, opt]=isp_makesonglist(varargin)
%ISP_MAKESONGLIST  Generate a song list from various collections
%
% SYNTAX
%   songArray=isp_makesonglist(dataset, options ...)
%   songArray=isp_makesonglist(options ...)
%
% DESCRIPTION
%   Reads track lists from either the ISMIR 2004 genre classification
%   training set, the ISMIR 2004 ballroom set, or from Dan Ellis'
%   artist20 or covers80 set, or a set of MIDI songs.
%
% INPUT
%   dataset:
%     Specifies the data set to read. Possible values are 'artist20',
%     'covers80', 'ismir2004ballroom', 'ismirgenre' and 'midi'.
%   options ...:
%     Structs or field/value pairs with the following field names:
%     artist20path, covers80path, ismir2004ballroompath, ismirgenrepath:
%       The root directory of the data collections. Defaults
%       are './artist20', './covers80', './' and './ismirgenre',
%       respectively.
%     midiInstruments:
%       For the midi set, specifies whether the same instrument are
%       used for all voices. Possible values are 'single' and 'multiple'.
%     midiset:
%       Use short (30 s) or long (600 s) midi files. Possible values:
%       'short', 'long'. Default: 'short'.
%     nInstruments:
%       For the midi set, specifies the maximum number of instrument
%       combinations.
%     nMelodies:
%       For the midi set, specifies the maximum number of melodies.
%
% HISTORY
%   Created by Jesper H. Jensen

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


    opt=struct('ismirgenrepath', './ismirgenre', ...
               'ismir2004ballroompath', './', ...
               'artist20path', './artist20', ...
               'covers80path', './covers80', ...
               'midiInstruments', 'single', ...
               'nInstruments', inf, ...
               'midiset', 'short', ...
               'nMidifiles', inf, ...
               'dataset', '');

    if ischar(varargin{1}) && ~isfield(opt, varargin{1})
        opt.dataset=varargin{1};
        opt=isp_interpretarguments(opt, varargin{2:end});
    else
        opt=isp_interpretarguments(opt, varargin{:});
    end        

    switch opt.dataset
      case 'midi'
        [songArray, opt]=midicollection(opt);
      case 'artist20'
        songArray=readArtist20tracklist(opt.artist20path);
      case 'covers80'
        songArray=readCoversongTracklist(opt.covers80path);
      case 'ismirgenre'
        songArray=readIsmirTracklist(opt.ismirgenrepath);
      case 'ismir2004ballroom'
        songArray=readIsmirBallroomTracklist(opt.ismir2004ballroompath);
      otherwise
        error(['Invalid data set ' dataset '.'])
    end
end



% readArtist20tracklist  Read song information for artist20 data set
%
% Syntax: songArray = readCoversongTracklist(coversongPath)
%
% Input:
%     coversongPath: Path to Dan Ellis' cover song set
%
% Output:
%     songArray: Array with information about the songs

function songArray = readArtist20tracklist(artist20path)

    songArray = [];
    for listNo = 0:5

        listFileName = fullfile(artist20path, 'lists', ...
                                ['a20-cut' num2str(listNo) '-tracks.list']);
        
        if ~exist(listFileName, 'file')
            error(['File ' listFileName ' does not exist'])
        end
        
        fid = fopen(listFileName, 'r');
        textCell=textscan(fid, '%s');
        fclose(fid);
        
        nSongs=length(textCell{1});
        for n=1:nSongs
            filename = textCell{1}{n};
            [artist, rest] = strtok([filename '//-'], '/');
            [album, rest] = strtok(rest(2:end), '/');
            [track, rest] = strtok(rest(2:end), '-');
            [title, rest] = strtok(rest(2:end), '/');

            song.filename = fullfile(artist20path, 'mp3s-32k', ...
                                     [filename '.mp3']);
            song.artist = artist;
            song.album = album;
            song.track = str2num(track);
            song.title = title;
            song.dataset = listNo;
            
            if ~exist(song.filename, 'file')
                error(['File ' song.filename ' not found.'])
            end

            if isempty(songArray)
                songArray = song;
            else
                songArray(end+1,1) = song;
            end
        end
    end
end




% readCoversongTracklist  Read song information for coversong data sets
%
% Syntax: songArray = readCoversongTracklist(coversongPath)
%
% Input:
%     coversongPath: Path to Dan Ellis' cover song set
%
% Output:
%     songArray: Array with information about the songs

function songArray = readCoversongTracklist(coverPath)

    songArray = [];
    for listNo = 1:2
        
        listFileName = fullfile(coverPath, 'covers32k', ['list' num2str(listNo) '.list']);
        
        if ~exist(listFileName, 'file')
            error(['File ' listFileName ' does not exist'])
        end
        
        fid = fopen(listFileName, 'r');
        textCell=textscan(fid, '%s');
        fclose(fid);
        
        nSongs=length(textCell{1});
        for n=1:nSongs
            filename = textCell{1}{n};
            [title, rest] = strtok([filename '/++-'], '/');
            [artist, rest] = strtok(rest(2:end), '+');
            [album, rest] = strtok(rest(2:end), '+');
            [track, rest] = strtok(rest(2:end), '-');
            [title2, rest] = strtok(rest(2:end), '/');

            %if ~isequal(title, title2)
            %    %error(['Error interpreting song name ' filename])
            %end

            song.filename = fullfile(coverPath, 'covers32k', [filename '.mp3']);
            song.artist = artist;
            song.album = album;
            song.track = str2num(track);
            song.title = title;
            song.dataset = listNo;
            
            if ~exist(song.filename, 'file')
                error(['File ' song.filename ' not found.'])
            end

            if isempty(songArray)
                songArray = song;
            else
                songArray(end+1,1) = song;
            end
        end
    end
end


function songArray = readIsmirBallroomTracklist(directory)
    datadir = fullfile(directory,'BallroomData', '');
    annodir=fullfile(directory,'BallroomAnnotations','ballroomGroundTruth','');
    filename = fullfile(datadir, 'allBallroomFiles');
    if ~exist(filename, 'file')
        error(['Cannot find file ' filename '.'])
    end
    
    fid=fopen(filename, 'r');

    filenames=textscan(fid, './%s\n');
    filenames = filenames{1};

    songArray(1:numel(filenames),1) = struct;
    for n=1:numel(filenames)
        songArray(n).filename = fullfile(datadir, filenames{n});
        if ~exist(songArray(n).filename, 'file')
            error(['Cannot find file ' filename '.'])
        end
        [tmp1, tmp2, tmp3] = fileparts(filenames{n});
        % The regexprep is to change Rumba-??? into Rumba
        songArray(n).style = regexprep(tmp1, '-.*', '');
        songArray(n).title = tmp2;
        gtfile = fullfile(annodir, [tmp2 '.bpm']);
        if ~exist(gtfile, 'file')
            error(['Cannot find file ' filename '.'])
        end
        gtfid = fopen(gtfile, 'r');
        songArray(n).tempo = fscanf(gtfid, '%f');
        fclose(gtfid);
    end


    
    %The test data definition is from http://mtg.upf.edu/ismir2004/contest/rhythmContest/TestData.txt
    testData = {'Albums-Cafe_Paradiso-06', 1
                'Albums-I_Like_It2-01', 1
                'Albums-Latin_Jam2-04', 1
                'Albums-Latin_Jam3-02', 1
                'Albums-Latin_Jam4-01', 1
                'Albums-Latin_Jam4-11', 1
                'Albums-Latin_Jam5-07', 1
                'Albums-Latino_Latino-04', 1
                'Albums-Macumba-01', 1
                'Albums-Mambo_Kings-03', 1
                'Albums-Mambo_Kings-10', 1
                'Albums-Pais_Tropical-06', 1
                'Media-103402', 1
                'Media-103407', 1
                'Media-103408', 1
                'Media-103412', 1
                'Media-103416', 1
                'Media-103417', 1
                'Media-103508', 1
                'Media-103603', 1
                'Media-103617', 1
                'Media-103716', 1
                'Media-103801', 1
                'Media-103908', 1
                'Media-104010', 1
                'Media-104106', 1
                'Media-104108', 1
                'Media-105417', 1
                'Media-105602', 1
                'Media-105603', 1
                'Media-106102', 1
                'Media-106104', 1
                'Media-106105', 1
                'Albums-Cafe_Paradiso-14', 2
                'Albums-Cafe_Paradiso-16', 2
                'Albums-Commitments-08', 2
                'Albums-Latin_Jam3-12', 2
                'Albums-Macumba-14', 2
                'Albums-Pais_Tropical-16', 2
                'Media-103515', 2
                'Media-103611', 2
                'Media-103612', 2
                'Media-103714', 2
                'Media-103813', 2
                'Media-103918', 2
                'Media-104017', 2
                'Media-104116', 2
                'Media-104118', 2
                'Media-105420', 2
                'Media-106115', 2
                'Media-106116', 2
                'Albums-AnaBelen_Veneo-11', 3
                'Albums-Ballroom_Classics4-18', 3
                'Albums-Ballroom_Classics4-19', 3
                'Albums-Chrisanne1-14', 3
                'Albums-Step_By_Step-15', 3
                'Albums-Step_By_Step-16', 3
                'Media-100615', 3
                'Media-103314', 3
                'Media-104216', 3
                'Media-104314', 3
                'Media-104317', 3
                'Media-104415', 3
                'Media-104417', 3
                'Media-104418', 3
                'Media-104716', 3
                'Media-104815', 3
                'Media-104915', 3
                'Media-104916', 3
                'Media-104918', 3
                'Media-105017', 3
                'Media-105120', 3
                'Media-105122', 3
                'Media-105517', 3
                'Media-105819', 3
                'Media-105911', 3
                'Albums-AnaBelen_Veneo-13', 4
                'Albums-Cafe_Paradiso-10', 4
                'Albums-Cafe_Paradiso-11', 4
                'Albums-Fire-10', 4
                'Albums-Latin_Jam-05', 4
                'Albums-Latin_Jam2-08', 4
                'Albums-Latin_Jam3-09', 4
                'Albums-Latin_Jam4-03', 4
                'Albums-Latin_Jam4-08', 4
                'Albums-Latin_Jam4-13', 4
                'Albums-Latino_Latino-05', 4
                'Albums-Macumba-09', 4
                'Albums-Pais_Tropical-12', 4
                'Media-103511', 4
                'Media-103512', 4
                'Media-103513', 4
                'Media-103608', 4
                'Media-103709', 4
                'Media-103811', 4
                'Media-103914', 4
                'Media-104005', 4
                'Media-104006', 4
                'Media-104112', 4
                'Media-104114', 4
                'Media-104115', 4
                'Media-105611', 4
                'Media-105612', 4
                'Media-105613', 4
                'Media-106110', 4
                'Albums-Cafe_Paradiso-02', 5
                'Albums-Cafe_Paradiso-03', 5
                'Albums-Fire-15', 5
                'Albums-Latin_Jam-06', 5
                'Albums-Latin_Jam2-03', 5
                'Albums-Latin_Jam2-14', 5
                'Albums-Latin_Jam3-04', 5
                'Albums-Latin_Jam3-06', 5
                'Albums-Latin_Jam4-02', 5
                'Albums-Latin_Jam4-12', 5
                'Albums-Latin_Jam5-01', 5
                'Albums-Pais_Tropical-01', 5
                'Albums-Pais_Tropical-03', 5
                'Media-103502', 5
                'Media-103505', 5
                'Media-103605', 5
                'Media-103705', 5
                'Media-103815', 5
                'Media-103902', 5
                'Media-105401', 5
                'Media-105416', 5
                'Media-105607', 5
                'Media-105608', 5
                'Media-106002', 5
                'Media-106004', 5
                'Media-106118', 5
                'Albums-Chrisanne2-02', 6
                'Albums-Commitments-11', 6
                'Albums-Secret_Garden-06', 6
                'Media-103302', 6
                'Media-103303', 6
                'Media-104204', 6
                'Media-104205', 6
                'Media-104301', 6
                'Media-104303', 6
                'Media-104402', 6
                'Media-104403', 6
                'Media-104503', 6
                'Media-104601', 6
                'Media-104703', 6
                'Media-104704', 6
                'Media-104801', 6
                'Media-104803', 6
                'Media-104902', 6
                'Media-104903', 6
                'Media-104904', 6
                'Media-105002', 6
                'Media-105003', 6
                'Media-105006', 6
                'Media-105101', 6
                'Media-105102', 6
                'Media-105103', 6
                'Media-105105', 6
                'Media-105213', 6
                'Media-105301', 6
                'Media-105304', 6
                'Media-105701', 6
                'Media-105801', 6
                'Media-105901', 6
                'Albums-Ballroom_Classics4-06', 7
                'Albums-Chrisanne2-04', 7
                'Albums-Chrisanne2-05', 7
                'Albums-Chrisanne2-06', 7
                'Albums-Chrisanne3-04', 7
                'Albums-Step_By_Step-05', 7
                'Albums-Step_By_Step-08', 7
                'Albums-StrictlyDancing_Tango-02', 7
                'Albums-StrictlyDancing_Tango-09', 7
                'Media-100605', 7
                'Media-100606', 7
                'Media-100607', 7
                'Media-103304', 7
                'Media-103305', 7
                'Media-104206', 7
                'Media-104707', 7
                'Media-104708', 7
                'Media-104709', 7
                'Media-104808', 7
                'Media-105109', 7
                'Media-105305', 7
                'Media-105505', 7
                'Media-105506', 7
                'Media-105508', 7
                'Media-105705', 7
                'Media-105904', 7
                'Albums-Ballroom_Classics4-11', 8
                'Albums-Ballroom_Classics4-12', 8
                'Albums-Ballroom_Classics4-14', 8
                'Albums-Chrisanne3-09', 8
                'Media-103307', 8
                'Media-104210', 8
                'Media-104212', 8
                'Media-104307', 8
                'Media-104309', 8
                'Media-104607', 8
                'Media-104809', 8
                'Media-104811', 8
                'Media-105010', 8
                'Media-105111', 8
                'Media-105307', 8
                'Media-105318', 8
                'Media-105509', 8
                'Media-105709', 8
                'Media-105809', 8
                'Media-105907', 8};

    % Training data definition is from
    % http://mtg.upf.edu/ismir2004/contest/rhythmContest/TrainingDataURLs.txt
    trainingData = {
        'Albums-Cafe_Paradiso-05', 1
        'Albums-Cafe_Paradiso-07', 1
        'Albums-Cafe_Paradiso-08', 1
        'Albums-Fire-01', 1
        'Albums-Fire-08', 1
        'Albums-Fire-14', 1
        'Albums-I_Like_It2-02', 1
        'Albums-Latin_Jam-01', 1
        'Albums-Latin_Jam-04', 1
        'Albums-Latin_Jam-07', 1
        'Albums-Latin_Jam-13', 1
        'Albums-Latin_Jam2-05', 1
        'Albums-Latin_Jam2-06', 1
        'Albums-Latin_Jam2-13', 1
        'Albums-Latin_Jam3-01', 1
        'Albums-Latin_Jam3-03', 1
        'Albums-Latin_Jam4-06', 1
        'Albums-Latin_Jam5-02', 1
        'Albums-Latin_Jam5-11', 1
        'Albums-Latino_Latino-01', 1
        'Albums-Latino_Latino-06', 1
        'Albums-Macumba-02', 1
        'Albums-Macumba-03', 1
        'Albums-Macumba-04', 1
        'Albums-Mambo_Kings-11', 1
        'Albums-Mambo_Kings-12', 1
        'Albums-Pais_Tropical-05', 1
        'Albums-Pais_Tropical-07', 1
        'Albums-Pais_Tropical-08', 1
        'Media-103401', 1
        'Media-103403', 1
        'Media-103404', 1
        'Media-103405', 1
        'Media-103406', 1
        'Media-103409', 1
        'Media-103410', 1
        'Media-103411', 1
        'Media-103413', 1
        'Media-103414', 1
        'Media-103415', 1
        'Media-103507', 1
        'Media-103509', 1
        'Media-103510', 1
        'Media-103601', 1
        'Media-103602', 1
        'Media-103604', 1
        'Media-103701', 1
        'Media-103702', 1
        'Media-103703', 1
        'Media-103704', 1
        'Media-103802', 1
        'Media-103803', 1
        'Media-103804', 1
        'Media-103816', 1
        'Media-103906', 1
        'Media-103907', 1
        'Media-103909', 1
        'Media-103910', 1
        'Media-104009', 1
        'Media-104011', 1
        'Media-104012', 1
        'Media-104107', 1
        'Media-104109', 1
        'Media-104110', 1
        'Media-104111', 1
        'Media-105405', 1
        'Media-105406', 1
        'Media-105407', 1
        'Media-105601', 1
        'Media-105604', 1
        'Media-105605', 1
        'Media-106005', 1
        'Media-106006', 1
        'Media-106007', 1
        'Media-106008', 1
        'Media-106101', 1
        'Media-106103', 1
        'Media-106117', 1
        'Albums-Cafe_Paradiso-15', 2
        'Albums-Commitments-10', 2
        'Albums-Fire-07', 2
        'Albums-Fire-12', 2
        'Albums-Latin_Jam-11', 2
        'Albums-Latin_Jam-12', 2
        'Albums-Latin_Jam2-11', 2
        'Albums-Latin_Jam2-12', 2
        'Albums-Latin_Jam3-11', 2
        'Albums-Latin_Jam4-05', 2
        'Albums-Latin_Jam4-10', 2
        'Albums-Latin_Jam5-04', 2
        'Albums-Latin_Jam5-09', 2
        'Albums-Macumba-15', 2
        'Albums-Macumba-16', 2
        'Albums-Pais_Tropical-13', 2
        'Albums-Pais_Tropical-14', 2
        'Albums-Pais_Tropical-15', 2
        'Media-103516', 2
        'Media-103517', 2
        'Media-103613', 2
        'Media-103614', 2
        'Media-103619', 2
        'Media-103713', 2
        'Media-103718', 2
        'Media-103719', 2
        'Media-103814', 2
        'Media-103819', 2
        'Media-103917', 2
        'Media-103919', 2
        'Media-103920', 2
        'Media-104015', 2
        'Media-104016', 2
        'Media-104117', 2
        'Media-104119', 2
        'Media-105413', 2
        'Media-105414', 2
        'Media-105415', 2
        'Media-106015', 2
        'Media-106016', 2
        'Media-106017', 2
        'Media-106114', 2
        'Albums-Ballroom_Classics4-20', 3
        'Albums-Ballroom_Magic-15', 3
        'Albums-Ballroom_Magic-16', 3
        'Albums-Ballroom_Magic-17', 3
        'Albums-Chrisanne1-12', 3
        'Albums-Chrisanne1-13', 3
        'Albums-Chrisanne2-11', 3
        'Albums-Chrisanne2-12', 3
        'Albums-Chrisanne3-13', 3
        'Albums-Chrisanne3-14', 3
        'Albums-Chrisanne3-15', 3
        'Albums-Step_By_Step-17', 3
        'Albums-Step_By_Step-18', 3
        'Media-100614', 3
        'Media-100616', 3
        'Media-103312', 3
        'Media-103313', 3
        'Media-104217', 3
        'Media-104218', 3
        'Media-104219', 3
        'Media-104315', 3
        'Media-104316', 3
        'Media-104416', 3
        'Media-104514', 3
        'Media-104515', 3
        'Media-104516', 3
        'Media-104615', 3
        'Media-104616', 3
        'Media-104617', 3
        'Media-104618', 3
        'Media-104717', 3
        'Media-104718', 3
        'Media-104816', 3
        'Media-104917', 3
        'Media-105016', 3
        'Media-105018', 3
        'Media-105019', 3
        'Media-105020', 3
        'Media-105119', 3
        'Media-105121', 3
        'Media-105207', 3
        'Media-105313', 3
        'Media-105314', 3
        'Media-105315', 3
        'Media-105320', 3
        'Media-105516', 3
        'Media-105518', 3
        'Media-105715', 3
        'Media-105716', 3
        'Media-105717', 3
        'Media-105816', 3
        'Media-105817', 3
        'Media-105818', 3
        'Media-105820', 3
        'Media-105912', 3
        'Media-105913', 3
        'Media-105914', 3
        'Albums-AnaBelen_Veneo-01', 4
        'Albums-AnaBelen_Veneo-03', 4
        'Albums-AnaBelen_Veneo-15', 4
        'Albums-Cafe_Paradiso-09', 4
        'Albums-Cafe_Paradiso-12', 4
        'Albums-Fire-03', 4
        'Albums-GloriaEstefan_MiTierra-01', 4
        'Albums-GloriaEstefan_MiTierra-04', 4
        'Albums-GloriaEstefan_MiTierra-06', 4
        'Albums-GloriaEstefan_MiTierra-08', 4
        'Albums-GloriaEstefan_MiTierra-11', 4
        'Albums-I_Like_It2-09', 4
        'Albums-Latin_Jam-02', 4
        'Albums-Latin_Jam-08', 4
        'Albums-Latin_Jam-14', 4
        'Albums-Latin_Jam2-07', 4
        'Albums-Latin_Jam2-09', 4
        'Albums-Latin_Jam2-15', 4
        'Albums-Latin_Jam3-07', 4
        'Albums-Latin_Jam3-08', 4
        'Albums-Latin_Jam5-03', 4
        'Albums-Latin_Jam5-08', 4
        'Albums-Latin_Jam5-12', 4
        'Albums-Latino_Latino-09', 4
        'Albums-Latino_Latino-10', 4
        'Albums-Macumba-10', 4
        'Albums-Macumba-11', 4
        'Albums-Macumba-12', 4
        'Albums-Pais_Tropical-09', 4
        'Albums-Pais_Tropical-10', 4
        'Albums-Pais_Tropical-11', 4
        'Media-103514', 4
        'Media-103609', 4
        'Media-103610', 4
        'Media-103615', 4
        'Media-103616', 4
        'Media-103710', 4
        'Media-103711', 4
        'Media-103712', 4
        'Media-103717', 4
        'Media-103809', 4
        'Media-103810', 4
        'Media-103812', 4
        'Media-103817', 4
        'Media-103818', 4
        'Media-103911', 4
        'Media-103912', 4
        'Media-103913', 4
        'Media-104007', 4
        'Media-104008', 4
        'Media-104113', 4
        'Media-105209', 4
        'Media-105408', 4
        'Media-105409', 4
        'Media-105410', 4
        'Media-105411', 4
        'Media-105412', 4
        'Media-105418', 4
        'Media-105419', 4
        'Media-105614', 4
        'Media-105615', 4
        'Media-106009', 4
        'Media-106010', 4
        'Media-106011', 4
        'Media-106012', 4
        'Media-106111', 4
        'Media-106112', 4
        'Media-106113', 4
        'Media-106119', 4
        'Albums-AnaBelen_Veneo-02', 5
        'Albums-Cafe_Paradiso-01', 5
        'Albums-Cafe_Paradiso-04', 5
        'Albums-Fire-02', 5
        'Albums-Fire-09', 5
        'Albums-Latin_Jam-03', 5
        'Albums-Latin_Jam-09', 5
        'Albums-Latin_Jam-15', 5
        'Albums-Latin_Jam2-01', 5
        'Albums-Latin_Jam2-02', 5
        'Albums-Latin_Jam3-05', 5
        'Albums-Latin_Jam4-07', 5
        'Albums-Latin_Jam5-06', 5
        'Albums-Latin_Jam5-10', 5
        'Albums-Latino_Latino-03', 5
        'Albums-Macumba-05', 5
        'Albums-Macumba-06', 5
        'Albums-Macumba-07', 5
        'Albums-Macumba-08', 5
        'Albums-Pais_Tropical-02', 5
        'Albums-Pais_Tropical-04', 5
        'Media-103501', 5
        'Media-103503', 5
        'Media-103504', 5
        'Media-103506', 5
        'Media-103606', 5
        'Media-103607', 5
        'Media-103618', 5
        'Media-103706', 5
        'Media-103707', 5
        'Media-103708', 5
        'Media-103715', 5
        'Media-103805', 5
        'Media-103806', 5
        'Media-103807', 5
        'Media-103808', 5
        'Media-103901', 5
        'Media-103903', 5
        'Media-103904', 5
        'Media-103905', 5
        'Media-104001', 5
        'Media-104002', 5
        'Media-104003', 5
        'Media-104004', 5
        'Media-104101', 5
        'Media-104102', 5
        'Media-104103', 5
        'Media-104104', 5
        'Media-104105', 5
        'Media-105402', 5
        'Media-105403', 5
        'Media-105404', 5
        'Media-105606', 5
        'Media-105609', 5
        'Media-106001', 5
        'Media-106003', 5
        'Media-106106', 5
        'Media-106107', 5
        'Media-106108', 5
        'Media-106109', 5
        'Albums-Ballroom_Classics4-01', 6
        'Albums-Ballroom_Classics4-02', 6
        'Albums-Ballroom_Classics4-03', 6
        'Albums-Ballroom_Classics4-04', 6
        'Albums-Ballroom_Classics4-05', 6
        'Albums-Ballroom_Magic-01', 6
        'Albums-Ballroom_Magic-02', 6
        'Albums-Ballroom_Magic-03', 6
        'Albums-Ballroom_Magic-04', 6
        'Albums-Ballroom_Magic-05', 6
        'Albums-Ballroom_Magic-18', 6
        'Albums-Chrisanne1-01', 6
        'Albums-Chrisanne1-02', 6
        'Albums-Chrisanne1-03', 6
        'Albums-Chrisanne2-01', 6
        'Albums-Chrisanne2-03', 6
        'Albums-Chrisanne3-01', 6
        'Albums-Chrisanne3-02', 6
        'Albums-Chrisanne3-03', 6
        'Albums-Fire-13', 6
        'Albums-Secret_Garden-01', 6
        'Albums-Secret_Garden-02', 6
        'Albums-Secret_Garden-05', 6
        'Albums-Step_By_Step-01', 6
        'Albums-Step_By_Step-02', 6
        'Albums-Step_By_Step-03', 6
        'Albums-Step_By_Step-04', 6
        'Media-100601', 6
        'Media-100602', 6
        'Media-100603', 6
        'Media-100604', 6
        'Media-103301', 6
        'Media-104201', 6
        'Media-104202', 6
        'Media-104203', 6
        'Media-104302', 6
        'Media-104401', 6
        'Media-104501', 6
        'Media-104502', 6
        'Media-104504', 6
        'Media-104602', 6
        'Media-104603', 6
        'Media-104604', 6
        'Media-104701', 6
        'Media-104702', 6
        'Media-104705', 6
        'Media-104802', 6
        'Media-104804', 6
        'Media-104805', 6
        'Media-104901', 6
        'Media-104905', 6
        'Media-105001', 6
        'Media-105004', 6
        'Media-105005', 6
        'Media-105007', 6
        'Media-105104', 6
        'Media-105106', 6
        'Media-105210', 6
        'Media-105211', 6
        'Media-105212', 6
        'Media-105214', 6
        'Media-105215', 6
        'Media-105302', 6
        'Media-105303', 6
        'Media-105316', 6
        'Media-105501', 6
        'Media-105502', 6
        'Media-105503', 6
        'Media-105504', 6
        'Media-105702', 6
        'Media-105703', 6
        'Media-105704', 6
        'Media-105802', 6
        'Media-105803', 6
        'Media-105804', 6
        'Media-105805', 6
        'Media-105902', 6
        'Albums-Ballroom_Classics4-07', 7
        'Albums-Ballroom_Classics4-08', 7
        'Albums-Ballroom_Classics4-09', 7
        'Albums-Ballroom_Classics4-10', 7
        'Albums-Ballroom_Magic-06', 7
        'Albums-Ballroom_Magic-07', 7
        'Albums-Ballroom_Magic-08', 7
        'Albums-Chrisanne1-04', 7
        'Albums-Chrisanne1-05', 7
        'Albums-Chrisanne1-06', 7
        'Albums-Chrisanne3-05', 7
        'Albums-Chrisanne3-06', 7
        'Albums-Step_By_Step-06', 7
        'Albums-Step_By_Step-07', 7
        'Albums-StrictlyDancing_Tango-01', 7
        'Albums-StrictlyDancing_Tango-03', 7
        'Albums-StrictlyDancing_Tango-04', 7
        'Albums-StrictlyDancing_Tango-05', 7
        'Albums-StrictlyDancing_Tango-06', 7
        'Albums-StrictlyDancing_Tango-07', 7
        'Albums-StrictlyDancing_Tango-08', 7
        'Albums-StrictlyDancing_Tango-10', 7
        'Albums-StrictlyDancing_Tango-11', 7
        'Albums-StrictlyDancing_Tango-12', 7
        'Albums-StrictlyDancing_Tango-13', 7
        'Albums-StrictlyDancing_Tango-14', 7
        'Albums-StrictlyDancing_Tango-15', 7
        'Media-103306', 7
        'Media-104207', 7
        'Media-104208', 7
        'Media-104304', 7
        'Media-104305', 7
        'Media-104306', 7
        'Media-104404', 7
        'Media-104405', 7
        'Media-104505', 7
        'Media-104506', 7
        'Media-104605', 7
        'Media-104606', 7
        'Media-104706', 7
        'Media-104710', 7
        'Media-104806', 7
        'Media-104807', 7
        'Media-104906', 7
        'Media-104907', 7
        'Media-105008', 7
        'Media-105009', 7
        'Media-105107', 7
        'Media-105108', 7
        'Media-105306', 7
        'Media-105317', 7
        'Media-105507', 7
        'Media-105706', 7
        'Media-105707', 7
        'Media-105806', 7
        'Media-105807', 7
        'Media-105808', 7
        'Media-105903', 7
        'Media-105905', 7
        'Media-105906', 7
        'Albums-Ballroom_Classics4-13', 8
        'Albums-Ballroom_Magic-09', 8
        'Albums-Ballroom_Magic-10', 8
        'Albums-Ballroom_Magic-11', 8
        'Albums-Chrisanne1-07', 8
        'Albums-Chrisanne1-08', 8
        'Albums-Chrisanne2-07', 8
        'Albums-Chrisanne3-07', 8
        'Albums-Chrisanne3-08', 8
        'Albums-Step_By_Step-09', 8
        'Albums-Step_By_Step-10', 8
        'Albums-Step_By_Step-11', 8
        'Media-100608', 8
        'Media-100609', 8
        'Media-100610', 8
        'Media-103308', 8
        'Media-103315', 8
        'Media-104209', 8
        'Media-104211', 8
        'Media-104308', 8
        'Media-104406', 8
        'Media-104407', 8
        'Media-104408', 8
        'Media-104409', 8
        'Media-104507', 8
        'Media-104508', 8
        'Media-104608', 8
        'Media-104609', 8
        'Media-104711', 8
        'Media-104712', 8
        'Media-104810', 8
        'Media-104908', 8
        'Media-104909', 8
        'Media-105011', 8
        'Media-105110', 8
        'Media-105112', 8
        'Media-105308', 8
        'Media-105309', 8
        'Media-105510', 8
        'Media-105511', 8
        'Media-105708', 8
        'Media-105710', 8
        'Media-105810', 8
        'Media-105811', 8
        'Media-105812', 8 };

    [titles{1:length(songArray),1}]=songArray.title;
    testData(:,3) = {'test'};
    trainingData(:,3) = {'training'};
    
    dataDesc = [testData; trainingData];
    [dummy, descidx] = sort(dataDesc(:,1));
    [dummy, songidx] = sort(titles);

    if ~all(strcmp(titles(songidx), dataDesc(descidx,1)))
        error('Data set has changed since this m-file was written.')
    end
    [songArray(songidx).dataset] = dataDesc{descidx,3};
    %[songArray(songidx).number] = dataDesc{descidx,2};
end





% readIsmirTracklist reads the ISMIR tracklist.csv and creates an array of structs with song information
%
% Syntax: songArray = readIsmirTracklist(directory)
%
% Input:
%     directory: The root directory where song information should be read from.
%         If it is not specified, a sensible default for the KOM network is
%         chosen.
%
% Output:
%     songArray: Array with information about the songs

function songArray = readIsmirTracklist(directory)
    filename = fullfile(directory,'Genre_part1','tracklist.csv');
    directories{1} = fullfile(directory,'Genre_part1','');
    directories{2} = fullfile(directory,'Genre_part2','');

    if ~exist(filename, 'file')
        error(['File ' filename ' does not exist'])
    end
    
    fid=fopen(filename, 'r');
    textCell=textscan(fid, '%s%s%s%s%s%s', 'delimiter', ',');
    fclose(fid);
    nSongs=length(textCell{1});

    temp=regexprep(textCell{1}, '"', '');
    [songArray(1:nSongs,1).genre]=temp{:};

    temp=regexprep(textCell{2}, '"', '');
    [songArray(1:nSongs).artist]=temp{:};

    temp=regexprep(textCell{3}, '"', '');
    [songArray(1:nSongs).album]=temp{:};
    
    temp=regexprep(textCell{4}, '"', '');
    [songArray(1:nSongs).title]=temp{:};
    
    temp=regexprep(textCell{5}, '"', '');
    [songArray(1:nSongs).track]=temp{:};
    
    temp=regexprep(textCell{6}, '"', '');
    temp=regexprep(temp, '^\./', '');
    [songArray(1:nSongs).filename]=temp{:};

    for n=1:nSongs
        songArray(n).filename=findfilelocation(songArray(n).filename, directories);
    end
end

% Determine in which directory a file resides. Used by 'readIsmirTracklist'
function truefilename=findfilelocation(filename, directories)
    % Take care of inaccuracies in tracklist.csv
    for fname={filename, strrep(filename, '''', '_'), ...
               regexprep(filename, '^tracks/', '')}
        for iDir=directories(:)'
            truefilename=fullfile(iDir{1}, fname{1});
            if exist(truefilename, 'file')
                return
            end
        end
    end
    error(['Cannot find file ' filename '.'])
end


function [songArray,opts]=midicollection(opts)
%

    % Define names of MIDI files
    switch opts.midiset
      case 'long'
        mididir=fullfile(isp_toolboxpath, 'longmidifiles', '');
      case 'short'
        mididir=fullfile(isp_toolboxpath, 'shortmidifiles', '');
      otherwise
        error('Unknown midi set.')
    end

    midifilenames = { ...
        '50s Rock.mid', ...
        'Argentina.mid', ...
        'Bach Prelude.mid', ...
        'Backbeat.mid', ...
        'Big Band.mid', ...
        'Blue Groove.mid', ...
        'Boogie.mid', ...
        'Caribbean.mid', ...
        'Chill.mid', ...
        'Country.mid', ...
        'Debussy.mid', ...
        'Eastern Europe.mid', ...
        'Funkeasy.mid', ...
        'Hungarian Folk.mid', ...
        'Italy.mid', ...
        'Jazz.mid', ...
        'Latin Five.mid', ...
        'Metal.mid', ...
        'Mexico.mid', ...
        'New Orleans.mid', ...
        'Pop Mellow.mid', ...
        'Punkarama.mid', ...
        'Rock I.mid', ...
        'Shufflemix.mid', ...
        'Slow Dance.mid', ...
        'South Africa.mid', ...
        'Strummin.mid', ...
        'Strut.mid', ...
        'Texas Swing.mid', ...
        'Waltz.mid'};

    midifilenames=regexprep(midifilenames, '^(.*)$', ...
                            [mididir filesep '$1']);
    if opts.nMidifiles <= length(midifilenames)
        opts.midifiles=midifilenames(1:opts.nMidifiles);
    else
        opts.midifiles=midifilenames;
        if ~isinf(opts.nMidifiles)
            warning('More midifiles requested than available.')
        end
    end

    instr=defineinstruments(opts.midiInstruments, opts.nInstruments);
    opts.instruments=instr;

    nFiles=numel(opts.midifiles);
    nInstr=numel(instr);

    nSongs = nFiles*nInstr;
    songArray=repmat(struct, nSongs, 1);

    kSong = 0;
    for iFile=1:nFiles
        [dummy, bname] = fileparts(midifilenames{iFile});
        for jInstr=1:nInstr
            kSong = kSong + 1;
            songArray(kSong).filename = midifilenames{iFile};
            songArray(kSong).melody = bname;
            songArray(kSong).instrument = instr{jInstr};
            songArray(kSong).modification.percussion = false;
        end
    end
    
end

    
function midiinstrument = defineinstruments(type, nInstruments)

    switch type

      case 'single'
        midiinstrument = num2cell([1 11 14 15 20 23 25 37 41 ...
                            47 53 54 57 66 71 74 76 77 79 81 82 85 89 93 ...
                            94 97 105 110 113 115]);
        if nInstruments > length(midiinstrument) && ~isinf(nInstruments)
            error('Unsupported number of selected instruments.')
        elseif nInstruments < length(midiinstrument)
            warning('Fewer instruments selected than expected.')
            midiinstrument = midiinstrument(1:nInstruments);
        end

      case 'random'
        if nInstruments > 112
            instruments = randperm(128); % Ignore percussive instruments
        else
            instruments = randperm(112); % Ignore percussive instruments
        end
        instruments = sort(instruments(1:nInstruments));
        midiinstrument = num2cell(instruments);

      case 'frommidifiles'
        error('Sorry, this functionality is not currently implemented.')
        midiinstrument = cell(1,length(midifile));
        for n=1:length(midifile)
            fprintf('Reading file %s.\n', midifile{n})
            [nmat, instrmat, ctrlmat] = ...
                readmidi2(midifile{n});
            midiinstrument{n} = instrmat( instrmat(:,2)~=10, 3);
            instrlength(n) = length(midiinstrument{n});
        end
        if ~all(instrlength == instrlength(1))
            error(['MIDI files have different number of non-percussive ' ...
                   'channels'])
        end

      case 'multiple'
        instrumentSets = {[25 27 31], [1 13 19], [33 36 59]};

        nInstrumentsPerSong = length(instrumentSets);
        for n=1:nInstrumentsPerSong
            nInstruments(n) = length(instrumentSets{n});
        end
        nInstrumentCombinations = prod(nInstruments);

        instrumentMatrix = zeros(nInstrumentsPerSong, nInstrumentCombinations);

        nCombinationsPerSong = nInstrumentCombinations;
        for n=1:nInstrumentsPerSong
            instr = instrumentSets{n}(:)';
            nInstr = length(instr);
            nCombinationsPerSong = nCombinationsPerSong / nInstr;
            tempMatrix = repmat(instr, nCombinationsPerSong, ...
                   nInstrumentCombinations / nCombinationsPerSong / nInstr);
            instrumentMatrix(n, :) = tempMatrix(:).';
        end
        
        midiinstrument = num2cell(instrumentMatrix, 1);

      otherwise
        error('Unknown instrument definition type')
    end
end