%ISP_RUNJOB  Execute pending, distributed jobs.
%
% SYNTAX
%   isp_runjob(jobDirectory)
%
% DESCRIPTION
%   Execute pending, distributed jobs.
%
% INPUT
%   jobDirectory:
%     Directory that contains the job description files. Default:
%     './jobDirectory'.
%
% SEE ALSO
%   isp_distribute.
%
% HISTORY
%   Copyright 2005-2008 Jesper Højvang Jensen.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function isp_runjob(jobDirectory)

    if ~exist('jobDirectory')
        jobDirectory='distributedJobs';
    end

    [temp,hostname]=unix('hostname');
    if isempty(hostname) || temp~=0
        error('This machine does not have a hostname.')
    end
    hostname=sscanf(hostname,'%s'); % Remove trailing newline.

    % Set file and path names
    jobsDir = fullfile(jobDirectory,'pendingJobs','');
    resultsDir = fullfile(jobDirectory,'results','');
    locksDir = fullfile(jobDirectory,'locks','');
    finishedJobsDir = fullfile(jobDirectory,'finishedJobs','');
    lockDir = fullfile(locksDir, hostname, '');
    if ~exist(lockDir, 'dir')
        mkdir(lockDir);
    end

    allJobsDone = false;
    while ~allJobsDone
        t0 = clock;
        randnState = randn('state');
        randState = rand('twister');

        % Hope the flush command will alleviate access problems
        if isunix
            system(['fs flush "' fullfile(jobsDir,'') '" 2> /dev/null']);
        end
        files=dir(fullfile(jobsDir, 'job*.mat'));

        % Stop if no more files
        foundFile = false;
        for fileNum = 1:length(files)
            if ~isempty(regexp(files(fileNum).name, '^job[0-9]{6}.mat$'))
                foundFile = true;
                break
            end
        end
        if ~foundFile
            allJobsDone = true;
            break
        end

        % File names
        jobsFile = fullfile(jobsDir, files(fileNum).name);
        resultsFile = fullfile(resultsDir, strrep(files(fileNum).name, 'job', 'results'));

        % Skip if the lock file exists
        lockfile=fullfile(lockDir, files(fileNum).name);
        if exist(lockfile, 'file')
            fprintf(1, ['Lock file %s already exists in lock directory. ' ...
                        'Continuing ...\n'], lockfile)
            pause(rand)
            continue
	end

        % Lock the first .mat file
        [success,msg]=copyfile(jobsFile, lockDir);
        if ~success
            warning(['Error copying file ' jobsFile ' to directory ' lockDir '.']);
            %       %%% DEBUG BEGIN
            %       fid=fopen('distributelogfile','a');
            %       fprintf(fid, 'Error copying file %s to directory %s.\nError message: %s\n',jobsFile, lockDir, msg);
            %       fclose(fid);
            %       %%% DEBUG END
            continue
        end
        delete(jobsFile);
        if isunix
            system(['fs flush "' fullfile(jobsDir,'') '" 2> /dev/null']);
        end
        %     %%% DEBUG BEGIN
        %     if exist(jobsFile, 'file')
        %       fid=fopen('distributelogfile','a');
        %       fprintf(fid, 'Error deleting file %s.\n', jobsFile);
        %       fclose(fid);
        %     end
        %     %%% DEBUG END

        % Do simulation
        fprintf('Running file %s\n', lockfile);
        workingdir=pwd;
        currentpath=path;
        settings=load(lockfile);
        if isfield(settings, 'randnState')
            randnState=settings.randnState;
        end
        if isfield(settings, 'randState')
            randState=settings.randState;
        end

        cd(settings.cwd);
        if (~isdeployed)
            addpath(settings.currentPath);
        end
        output=evalfunction(settings);
        cd(workingdir);
        if (~isdeployed)
            path(currentpath);
        end
        save(resultsFile, 'output', 'randnState', 'randState')

        % Indicate job is finished
        [success, msg]=movefile(lockfile, finishedJobsDir);
        if success
            system(sprintf('echo %s finished %s in %d seconds >> "%s"', ...
                           hostname, files(fileNum).name, ...
                           etime(clock, t0), ...
                           fullfile(jobDirectory, 'hostsstatus')));
        else
            warning(['Error moving file ' lockfile ' to directory ' finishedJobsDir '.']);
            fprintf('Error message: %s\n', msg);
            system(sprintf('echo %s had error moving file %s to finished jobs directory. Error message: %s >> %s', ...
                           hostname, files(fileNum).name, msg, ...
                           fullfile(jobDirectory, 'hostsstatus')));
            %       %%% DEBUG BEGIN
            %       fid=fopen('distributelogfile','a');
            %       fprintf(fid, 'Error moving file %s to directory %s.\nError message: %s\n',lockfile ,finishedJobsDir , msg);
            %       fclose(fid);
            %       %%% DEBUG END
            continue
        end

    end
    fprintf('\nNo more pending jobs.\n')
end


function output=evalfunction(settings)

    % Set states of random number generators
    rand('twister', settings.randState);
    randn('state', settings.randnState);

    % Set global variables
    for globalVar=fieldnames(settings.globalVariables)'
        if strcmp(globalVar{1}, 'globalVar') || strcmp(globalVar{1}, 'settings') || strcmp(globalVar{1}, 'output')
            error(['Global variable ' globalVar{1} ' is causing a name conflict']);
        end
        eval(['global ' globalVar{1} '; ' globalVar{1} '= settings.globalVariables.' globalVar{1} ';' ]);
    end

    % Execute distributed command
    [output{1:settings.numLHSarguments}]= ...
        feval(settings.functionName, settings.inputArguments{:});
end
