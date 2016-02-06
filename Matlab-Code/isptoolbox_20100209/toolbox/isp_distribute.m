% ISP_DISTRIBUTE  Distribute job for parallel execution on multiple machines.
%
% SYNTAX
%   isp_distribute(command, jobDirectory)
%   isp_distribute(commandstring, nOutputs, functionNames)
%
% DESCRIPTION
%   Either schedule command for delayed and/or parallel execution, or
%   read the results of a previously scheduled jobs.
%
%   When first calling this function, you must specify the command 'start',
%   optionally followed by the directory which shall be used for temporary
%   files.
%
%   Next, call this function from a script/function with all the command
%   lines that are to be distributed to other computers (for instance a
%   test with different parameters). After this, either call
%   executependingjobs.m from a number of Matlab sessions, or use the
%   'distribute'-script to batch-process the jobs.
%
%   After completing all jobs, call this function with the command 'read'.
%   When the script/function is now run a second time (without changes!), the
%   results of the distributed commands will be loaded.
%
%   There are a few things to be aware of: The random number generators for
%   rand and randn are set to 'twister' and 'state', respectively, and their
%   states are saved when distributing a job. Each generator is advanced by
%   10000 numbers for each distributed task. If this is not enough to ensure
%   random experiments, they should manually be advanced between each call
%   to isp_distribute.  Each time 'distribute' is called with a job, it
%   saves all global variables and function arguments. A large global
%   variable or a large arguments might therefore take up a lot of harddisk
%   space.  Errors might occur if there are relative paths in the matlab
%   path. That is, instead of for instance addpath('utilities'), use
%   addpath(fullfile(pwd, 'utilities', '')).
%
% INPUT
%   command:
%     One of the following strings:
%     'start':
%       Specifies that the next jobs are to be distributed
%     'read':
%       Specifies that the results of the next jobs are to be read.
%     'execute':
%       Specifies that the next jobs are to be executed immediately.
%     'ask':
%       Ask whether to 'start', 'read' or 'execute'.
%     'skip':
%       Increment job number without actually reading results.
%     'state':
%       Change state number. Used to have parallel distributed job
%       sessions.
%   jobDirectory:
%     Directory that will contain the job description files, results
%     etc. Default: './jobDirectory'.
%   commandstring:
%     A Matlab command on one of the forms
%     'someFunction(arg1, arg2, ...)',
%     'variable1=someFunction',
%     '[variable1, variable2]=someFunction',
%     'variable1=someFunction(arg1, arg2, ...)', or
%     '[variable1, variable2]=someFunction(arg1, arg2, ...)'
%   nOutputs:
%     Number of output parameters. If this argument is not present or is the
%     empty matrix, the number of outputs are estimated. It is only
%     necessary to specify the number of outputs if the estimation
%     fails, which will happen if e.g. the "[outarg{1:N}] = somefunction"
%     syntax is used
%   functionNames: Cell array of function names that are called
%     indirectly by the function specified in 'commandstring' using
%     e.g. the feval function. Specifying this argument is only necessary
%     if the mcc compiler is used and it cannot automatically identify
%     all dependencies.
%
% OUTPUT
%   No arguments are returned directly. However, if the specified job has been
%   executed, the assignment specified in 'commandstring' will be
%   performed in the callers workspace.
%
% EXAMPLE
%   Create a script with the following contents:
%     for n=1:9
%       B=eye(8);
%       distribute('A{n}=pinv(B)');
%     end
%     whos
%
%   Run distributejob('start') and then run the script. This will create job
%   description files. Now call executependingjobs.m to execute the
%   jobs. Then run distributejob('read') and run the script again. The
%   pre-calculated results will now be loaded from disk.
%
% SEE ALSO
%   isp_runjob.
%
% HISTORY
%   Copyright 2005-2008 Jesper Højvang Jensen.
%   Last modified 01-01-2007.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


function isp_distribute(commandstring, nOutputs, functionNames)


    persistent nTimesCalled
    persistent specifyJobs
    persistent jobDirectory
    persistent activeState

    if isempty(activeState), activeState = 1; end

    % Handle input arguments
    switch commandstring
      case 'state'
        activeState = nOutputs;
        return
      case {'start', 'read', 'execute', 'ask'}
        if exist('nOutputs')
            jobDirectory{activeState}=nOutputs;
        elseif isempty(jobDirectory) || length(jobDirectory) < activeState ...
                || ~ischar(jobDirectory{activeState})
            jobDirectory{activeState}='distributedJobs';
        end
      otherwise
        if isempty(jobDirectory) || isempty(jobDirectory{activeState})
            error('Before specifying jobs or loading results, you must give the ''start'', ''read'', ''ask'' or ''execute'' parameter');
        end
    end

    if ~exist('functionNames')
        functionNames={};
    end

    if isstr(functionNames)
        functionNames={functionNames};
    end

    % Set file and path names
    jobsDir=fullfile(jobDirectory{activeState},'pendingJobs','');
    resultsDir=fullfile(jobDirectory{activeState},'results','');
    locksDir=fullfile(jobDirectory{activeState},'locks','');
    finishedJobsDir=fullfile(jobDirectory{activeState},'finishedJobs','');


    % Ask user what to do
    if strcmp(commandstring, 'ask')
        fprintf(1, '\n1. Start distributing jobs.\n')
        fprintf(1, '2. Read distributed jobs.\n')
        fprintf(1, '3. Execute jobs immediately.\n')
        answer = str2double(input('What do you want to do? ', 's'));
        switch(answer)
          case 1
            commandstring = 'start';
          case 2
            commandstring = 'read';
          case 3
            commandstring = 'execute';
          otherwise
            error('Invalid choice')
        end
    end


    % Interpret commands
    switch (commandstring)
      case 'start'
        nTimesCalled(activeState) = 0;
        specifyJobs{activeState} = 'specify';
        for dir={jobDirectory{activeState}, jobsDir, locksDir, finishedJobsDir, resultsDir}
            if ~exist(dir{1}, 'dir')
                mkdir(dir{1});
            end
        end
        % Set random number generators to the ones we query later on.
        %fprintf(['Setting random number generator for rand and randn to ' ...
        %         '''twister'' and ''state'', respectively.\n']);
        rand('twister', rand('twister'));
        randn('state', randn('state'));
        return
      case 'read'
        nTimesCalled(activeState) = 0;
        specifyJobs{activeState} = 'load';
        return
      case 'skip'
        nTimesCalled(activeState) = nTimesCalled(activeState) + 1;
        fprintf('Skipping job %s.\n', fullfile(jobsDir, sprintf('job%06d.mat', nTimesCalled(activeState))));
        return
      case 'execute'
        nTimesCalled(activeState) = 0;
        specifyJobs{activeState} = 'execute';
        return
      otherwise
    end

    % Interpret Matlab-command
    [functionName, arguments, lhs, numLHSarguments]=interpretCmdString(commandstring);
    if exist('nOutputs') && ~isempty(nOutputs)
        numLHSarguments=nOutputs;
    end


    nTimesCalled(activeState) = nTimesCalled(activeState) + 1;

    switch specifyJobs{activeState}
      case 'specify'
        % Create job file
        jobFile=fullfile(jobsDir, sprintf('job%06d.mat', nTimesCalled(activeState)));
        fprintf('Creating job file %s.\n', jobFile);
        globalVariables=findGlobalVariables;
        inputArguments={};
        for n=1:length(arguments)
            inputArguments{n}=evalin('caller',arguments{n});
        end
        cwd=pwd;
        currentPath=cleanpath;
        randnState=randn('state');
        randState=rand('twister');
        save(jobFile, 'cwd', 'globalVariables', 'functionName', 'numLHSarguments', 'inputArguments', 'currentPath', 'functionNames', 'randnState', 'randState');
        randn(100);
        rand(100);
      case 'load'
        % Load pre-calculated results
        resultsFile=fullfile(resultsDir, sprintf('results%06d.mat', nTimesCalled(activeState)));
        fprintf('Loading results from file %s.\n', resultsFile);
        load(resultsFile)
        if ~isempty(output)
            if ~isempty(whos('global','temporaryGlobalVariable')) || ...
                    evalin('caller','exist(''temporaryGlobalVariable'')')
                error('''temporaryGlobalVariable'' already exists. I need it.');
            end

            global temporaryGlobalVariable
            evalin('caller', 'global temporaryGlobalVariable');
            temporaryGlobalVariable=output;
            % Octave doesn''t support the 'global' option to clear
            % (version 2.9.14), so clearing 'temporaryGlobalVariable' in
            % the caller workspace is necessary.
            evalin('caller', [lhs '=deal(temporaryGlobalVariable{:}); ' ...
                                'clear temporaryGlobalVariable']);
            clear global temporaryGlobalVariable
        end
        randn(100);
        rand(100);
      case 'execute'
        evalin('caller', commandstring);
      otherwise
        error('This is not supposed to happen. There is a bug in the code.')
    end

end


function globalVariables=findGlobalVariables
    % Return a structure containing all global variables
    
    globalVariables=struct();
    for globalVariableName=who('global')';
        if strcmp(globalVariableName{1},'globalVariableNames') || ...
                strcmp(globalVariableName{1},'globalVariables')
            warning(['You have a global variable named ' globalVariableName{1} '. This variable name is also used for other purposes.'])
            continue
        end
        eval(['global ' globalVariableName{1}])
        globalVariables=setfield(globalVariables, globalVariableName{1}, eval(globalVariableName{1}));
    end
end


function [functionName, arguments, lhs, numLHSarguments]=interpretCmdString(commandstring)
    % Interpret command string

    assignmentPos=strfind(commandstring,'=');
    if length(assignmentPos)>1
        error('Command string contains more than one =.')
    end

    if length(assignmentPos)==0
        lhs='';
        rhs=strtrim(commandstring);
    else
        lhs=strtrim(commandstring(1:assignmentPos-1));
        rhs=strtrim(commandstring(assignmentPos+1:end));
    end

    argStartpos=strfind(rhs,'(');
    if isempty(argStartpos)
        functionName=rhs;
        arguments={};
    else
        functionName=strtrim(rhs(1:argStartpos(1)-1));
        argumentString=strtrim(rhs(argStartpos(1)+1:end));
        round=0;
        brackets=0;
        curly=0;
        argNum=1;
        arguments{argNum}='';
        for n=1:length(argumentString)
            switch argumentString(n)
                case '('
                    round=round+1;
                case ')'
                    round=round-1;
                case '['
                    brackets=brackets+1;
                case ']'
                    brackets=brackets-1;
                case '{'
                    curly=curly+1;
                case '}'
                    curly=curly-1;
            end

            if argumentString(n)==',' && round==0 && brackets==0 && curly==0
                argNum=argNum+1;
                arguments{argNum}='';
            else
                arguments{argNum}(end+1)=argumentString(n);
            end
        end
        if round~=-1 || curly~=0 || brackets~=0
            error('Error interpreting function arguments.');
        end
        rightparpos=strfind(arguments{argNum},')');
        arguments{argNum}=arguments{argNum}(1:rightparpos(end)-1);
    end

    if isempty(lhs)
        numLHSarguments=0;
    elseif isempty(strfind(lhs,'['))
        numLHSarguments=1;
    else
        brackPos=strfind(lhs,'[');
        insideBrackets=strtrim(lhs(brackPos(1)+1:end));
        insideBrackets=strrep(insideBrackets,',',' ');
        round=0;
        brackets=0;
        curly=0;
        argNum=1;
        for n=2:length(insideBrackets)
            switch insideBrackets(n)
                case '('
                    round=round+1;
                case ')'
                    round=round-1;
                case '['
                    brackets=brackets+1;
                case ']'
                    brackets=brackets-1;
                case '{'
                    curly=curly+1;
                case '}'
                    curly=curly-1;
            end

            if insideBrackets(n)==' ' && insideBrackets(n-1)~=' ' && round==0 && brackets==0 && curly==0
                argNum=argNum+1;
            end
        end
        if round~=0 || curly~=0 || brackets~=-1
            error('Error interpreting left hand side arguments.');
        end
        numLHSarguments=argNum;
    end
end


% Return the path with all entries starting with 'matlabroot' removed.
function cp=cleanpath
    if exist('matlabroot')
        persistent prevCp prevP
        p=path;
        if strcmp(p, prevP)
            cp = prevCp;
            return
        end
        prevP = p;
        cp='';
        root=matlabroot;
        while ~isempty(p)
            [t,p] = strtok(p, pathsep);
            if length(t) < length(root) || ~all(t(1:length(root))==root)
                cp=[cp pathsep t];
            end
        end
        if ~isempty(cp)
            cp(1)=[];
        end
        prevCp = cp;
    else
        cp=path;
    end
end