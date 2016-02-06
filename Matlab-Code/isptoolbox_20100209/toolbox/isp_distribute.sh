#!/bin/bash

echo Running $0

# Script to distribute matlab jobs to several computers.
#
# Copyright 2005-2006 Jesper Højvang Jensen
# Last modified 05-04-2006


# Use mcc to compile the matlab code? Using compiled code saves matlab
# licenses, but it might not work and it makes troubleshooting more difficult.
USEMCC=false
#USEMCC=true

# Specify if restrictions apply to the time interval where the
# distributed jobs are executed.

#RESTRICTTIME=true
RESTRICTTIME=false
STARTTIME=1000
STOPTIME=0745

# The following specifies if jobs shall terminate if another user
# starts a thread on the machine
STOPFOROTHERUSERS=true
#STOPFOROTHERUSERS=false

# Log screen output to locks/machine/logfile.$PID
LOG=true
#LOG=false

# The number of Matlab licenses to ensure always is left.
MINIMUMMATLABLICENSES=20

MATLABPATH="/afs/ies.auc.dk/user/jhj/no_backup/software/matlab-r2007b"
DISTRIBUTESCRIPTPATH=`dirname \`which $0\``
DEFAULTDISTRIBUTEDDIR='distributedJobs'



# Specifies whether we are running under sunOS or not
case `uname -m` in
    (i*86)
      SUNOS=false;;
    (sun*)
      SUNOS=true;;
esac

function validtime() {
# Test if we are still within then allowed time range
    if ! $RESTRICTTIME; then
      return;
    fi
    start=$STARTTIME
    stop=$STOPTIME
    curtime=`date +'%H%M'`
    if [ $((10#$start)) -gt $((10#$stop)) ]; then
      stop=$(( 10#$stop + 2400 ));
      if [ $((10#$curtime)) -lt $((10#$start)) ]; then
        curtime=$(( 10#$curtime + 2400 ))
      fi
    fi
#    echo times: $start $stop $curtime
    if [ $((10#$curtime)) -gt $((10#$start)) ] && [ $((10#$curtime)) -lt $((10#$stop)) ]; then
	true
    else
	false
    fi
}

function enoughlicensesleft() {
# Print how many matlab licenses are left
  totallicenses=`$MATLABPATH/etc/lmstat -i MATLAB|grep '^MATLAB [ ].*.[0-9]*. [ ]*'|tail -1|awk '{print $3}'`
  usedlicenses=`$MATLABPATH/etc/lmstat -f MATLAB|grep 'Users of MATLAB:  (Total of .* licenses issued;  Total of .* licenses in use)'|awk '{print $11}'`

  if  [ $(( totallicenses )) -gt 0 ] && [ $(( totallicenses - usedlicenses )) -lt $MINIMUMMATLABLICENSES ]; then
      false
  else
      true
  fi
}

function showthreads() {
# Print user threads. This is used to check if other users wish to run
# software on this computer
  echo `ps -eo user,comm|egrep -v "^ *(root|$USER|gdm|daemon|nobody|cupsys|postfix|syslog|news|man|hal|klog|identd) "`
}


function compile() {
    echo Compiling .m files

    # Clear the cache. Might otherwise cause trouble with the DISPLAY environment variable
    $MATLABPATH/bin/mcc -rmcache

    # It is necessary to do this cumbersome trick, because we want mcc to
    # inherit the path statement from matlab.
    tempdir="./distributedmatlabjobstemporarydirectory"
    DISPLAY= "$MATLABPATH/bin/matlab" -nodisplay -nojvm <<EOF
      addpath('$DISTRIBUTESCRIPTPATH');
      distdir='$distributedDir';
      jobdir=fullfile(distdir, 'pendingJobs', '');
      files=dir(fullfile(jobdir, 'job*.mat'));
      funcNames={};
      fprintf('Investigating job files\n');
      for fileNum=1:length(files)
          s=load(fullfile(jobdir, files(fileNum).name));
          if isempty(regexp(which(s.functionName), '.bi\$|^built-in'))
            funcNames{end+1}=s.functionName;
          end
      end
      funcNames=unique(funcNames);
      path(s.currentPath)
      fprintf('Compiling files\n');
      disp(funcNames)
%      mcc('-v','-m','-R','nojvm','-R','-nodisplay','-d',distdir,'executependingjobs',funcNames{:}, s.functionNames{:})
% Because of a bug in mbuild, we have to compile files to $tempdir and later move them to their right location
      mkdir('$tempdir')
      mcc('-m','-R','-nojvm','-R','-nodisplay','-d','$tempdir','isp_runjob',funcNames{:}, s.functionNames{:})
      fprintf('Done compiling\n');
      quit
EOF
      mkdir "$distributedDir/`uname -m`" >& /dev/null
      mv -f $tempdir/* "$distributedDir/`uname -m`"
      rmdir $tempdir
}

function startMatlab() {
    echo `hostname` StartMatlab: `pwd` $distributedDir
    if ! $USEMCC; then
	$MATLABPATH/bin/matlab -nojvm -nodisplay <<EOF
          addpath('$DISTRIBUTESCRIPTPATH');
          isp_runjob('$distributedDir')
          quit
EOF
    else
	if $SUNOS; then
	    # Sun OS
	    LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MATLABPATH/bin/sol2/"
	    export LD_LIBRARY_PATH
	    "$distributedDir/`uname -m`/isp_runjob" "$distributedDir"
	else
	    # Intel machine
	    LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MATLABPATH/bin/glnx86/:$MATLABPATH/sys/os/glnx86/:$distributedDir"
	    export LD_LIBRARY_PATH
	    LD_PRELOAD="`echo $MATLABPATH/sys/os/glnx86/libstdc++.so.?`" "$distributedDir/`uname -m`/isp_runjob" "$distributedDir"
	fi
    fi
}


function compiledistribute() {
# Distribute jobs to all 
    if $USEMCC; then
	compile
    fi

    distribute "$@"
}

function distribute() {
# Distribute jobs to all 
    if [ $# -lt 1 ]; then
	echo No hosts specified
	exit
    fi

    sleepTime=0
    for f in "$@"; do
# I don't remember why /afs/ies.auc.dk/sw/released/bin was explicitly added to PATH for the SUN servers. If it causes trouble not to do so, uncomment the following lines.
#	if $SUNOS; then
#	    echo sleep $sleepTime';' cd \"$PWD\"';' PATH='$PATH:/afs/ies.auc.dk/sw/released/bin' nice -19 screen -d -m \"$0\" internalfunction -d \"$distributedDir\" | ssh -x $f &
#	else
	echo sleep $sleepTime';' cd \"$PWD\"';' LOG=$LOG nice -19 screen -L -d -m bash -c \'sleep 5\; \"$0\" internalfunction -d \"$distributedDir\"\' | ssh -x $f &
	    # echo sleep $sleepTime';' cd \"$PWD\"';' nice -19 screen -d -m \"$0\" internalfunction -d \"$distributedDir\" '|' ssh -x $f >> testlog
#	fi
        sleepTime=$((sleepTime + 10))
    done
}

function internalfunction() {
    screenpid=$PPID
    lockpath="$distributedDir/locks/`hostname`"
    settingsfile="$lockpath/settings.cfg"
    logfile="$lockpath/screenlog.$screenpid"

    if ! [ -d $lockpath ]; then
	mkdir "$lockpath"
    fi

    # Store current settings
    echo "# Settings for jobs running on `hostname`" > "$settingsfile"
    for var in RESTRICTTIME STARTTIME STOPTIME STOPFOROTHERUSERS LOG MINIMUMMATLABLICENSES; do
	echo "$var=${!var}" >> "$settingsfile"
    done

    oldLogSetting=false

    startMatlab -d "$distributedDir" &
    matlabPID=$!

    echo `hostname` `date` PID $screenpid: Starting >> "$distributedDir/hostsstatus"

    jobs=`showthreads`
    while true; do
      # Read settings
      if [ -f "$settingsfile" ]; then
	  source "$settingsfile"
      fi

      # Update log setting if necessary
      if [ "$LOG" != "$oldLogSetting" ]; then
	  oldLogSetting="$LOG"
	  if $LOG; then
	      screen -X logfile "$logfile"
	      screen -X log on
	      screen -X logtstamp on
	      echo `hostname` `date` PID $screenpid: Logging screen output to $logfile >> "$distributedDir/hostsstatus"
	  else
	      screen -X log off
	      screen -X logtstamp off
	      echo `hostname` `date` PID $screenpid: Stopped logging screen output >> "$distributedDir/hostsstatus"
	  fi
      fi
      

      if $STOPFOROTHERUSERS && [ "$jobs" != "`showthreads`" ]; then
	echo n|mv "$lockpath"/job*.mat "$distributedDir/pendingJobs"
	echo `hostname` $screenpid killed because another user needed the machine >> "$distributedDir/hostsstatus"
	echo Old jobs: $jobs >> "$distributedDir/hostsstatus"
	echo New jobs: `showthreads` >> "$distributedDir/hostsstatus"
	kill -9 $screenpid
      fi
      if ! validtime ; then
	echo n|mv "$lockpath"/job*.mat "$distributedDir/pendingJobs"
	echo `hostname` $screenpid killed because dawn is breaking >> "$distributedDir/hostsstatus"
	kill -9 $screenpid
      fi
      if ! $USEMCC && ! enoughlicensesleft ; then
	echo n|mv "$lockpath"/job*.mat "$distributedDir/pendingJobs"
	echo `hostname` $screenpid killed because too many matlab licenses are in use >> "$distributedDir/hostsstatus"
	kill -9 $screenpid
      fi
      if [ -z "`ps -p $matlabPID|grep -v PID`" ]; then
	echo `hostname` $screenpid have finished >> "$distributedDir/hostsstatus"
	if $LOG; then
	    sleep 60 # Wait for screen to flush output to log file
	fi
	return 0
      fi
      sleep 30
    done

#    wait $screenpid
#    if [ $? -ne 0 ]; then
#      echo `hostname` $screenpid cancelled >> "$distributedDir/hostsstatus"
#    fi
      
#    echo `hostname` $screenpid finished >> "$distributedDir/hostsstatus"
}

function clear() {
    rm "$distributedDir/hostsstatus"
}

function killjobs() {
    for c in "$@"; do
        ssh -x $c 'pkill -9 distribute; sleep 1; screen -wipe' &
	lockpath="$distributedDir/locks/$c"
	echo n|mv "$lockpath"/job*.mat "$distributedDir/pendingJobs"
    done

}


if [ $# -eq 0 ]; then
    echo "Start a distributed matlab job on other machines"
    echo
    echo "Syntax: $0 command [options1] [option2] ..."
    echo
    echo "'command' must be one of the following:"
    echo
    echo "clear: clear list of used computers"
    echo "compile: Use mcc to compile standalone job"
    echo "killjobs: Kill all jobs on specified hosts. If no hosts are specified, kill all jobs."
    echo "run: Start matlab on this computer and begin executing jobs."
    echo "compiledistribute: Compile if necessary and distribute jobs to the computers specified by the following options."
    echo "distribute: Distributes the jobs to the computers specified by the following options without compiling."
    echo
    echo "Options:"
    echo " -d directory   The directory containing information about the matlab jobs to execute. It is also used for bookkeeping. This option must be specified right after the command"
    echo
    exit
fi

distributedDir="$DEFAULTDISTRIBUTEDDIR"

if [ $# -gt 2 ]; then
    if [ "$2" == '-d' ]; then
	command="$1"
	distributedDir="`echo "$3"|sed -e 'sa/$aa'`"
	shift 3
    elif [ "$1" == "-d" ]; then
	distributedDir="`echo "$2"|sed -e 'sa/$aa'`"
	command="$3"
	shift 3
    fi
else
    command="$1"
    shift
fi

if [ ! -d "$distributedDir" ]; then
    echo "Error: Directory $distributedDir does not exist."
    exit
fi

case $command in
    (compile)
        compile;;
    (run)
        startMatlab;;
    (startMatlab) # This is only called internally.
        startMatlab;;
    (compiledistribute)
        compiledistribute "$@";;
    (distribute)
        distribute "$@";;
    (killjobs)
        if [ $# -eq 0 ]; then 
	    killjobs `cat "$distributedDir/hostsstatus"|awk '{print $1}'|sort|uniq`
	else
	    killjobs "$@"
	fi ;;
    (internalfunction)
        internalfunction;;
    (clear)
        clear;;
    (*)
        echo "Invalid command specified";;
esac
