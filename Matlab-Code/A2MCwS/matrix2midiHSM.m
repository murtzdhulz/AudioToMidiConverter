function midi=matrix2midiHSM(M,ticks_per_quarter_note,timesig)

if nargin < 2
  ticks_per_quarter_note = 300;
end

if nargin < 3
  timesig = [4,2,24,8];
end

 midi.format = 0;
 midi.ticks_per_quarter_note = ticks_per_quarter_note;

tempo = 500000; 
trM = M;
note_events_onoff = [];
note_events_n = [];
note_events_ticktime = [];

for j=1:size(trM,1)
    % note on event:
    note_events_onoff(end+1)    = 1;
    note_events_n(end+1)        = j;
    note_events_ticktime(end+1) = 1e6 * trM(j,5) * ticks_per_quarter_note / tempo;
    
    % note off event:
    note_events_onoff(end+1)    = 0;
    note_events_n(end+1)        = j;
    note_events_ticktime(end+1) = 1e6 * trM(j,6) * ticks_per_quarter_note / tempo;
  end
 msgCtr = 1;
  
  % set tempo...
  midi.track(1).messages(msgCtr).deltatime = 0;
  midi.track(1).messages(msgCtr).type = 81;
  midi.track(1).messages(msgCtr).midimeta = 0;
  midi.track(1).messages(msgCtr).data = encode_int(tempo,3);
  midi.track(1).messages(msgCtr).chan = [];
  msgCtr = msgCtr + 1;
  
  % set time sig...
  midi.track(1).messages(msgCtr).deltatime = 0;
  midi.track(1).messages(msgCtr).type = 88;
  midi.track(1).messages(msgCtr).midimeta = 0;
  midi.track(1).messages(msgCtr).data = timesig(:);
  midi.track(1).messages(msgCtr).chan = [];
  msgCtr = msgCtr + 1;
  
  [junk,ord] = sort(note_events_ticktime);

  prevtick = 0;
  for j=1:length(ord)
    
    n = note_events_n(ord(j));
    cumticks = note_events_ticktime(ord(j));
    
    midi.track(1).messages(msgCtr).deltatime = cumticks - prevtick;
    midi.track(1).messages(msgCtr).midimeta = 1; 
    midi.track(1).messages(msgCtr).chan = trM(n,2);
    midi.track(1).messages(msgCtr).used_running_mode = 0;

    if (note_events_onoff(ord(j))==1)
      % note on:
      midi.track(1).messages(msgCtr).type = 144;
      midi.track(1).messages(msgCtr).data = [trM(n,3); trM(n,4)];
    else
      %-- note off msg:
      %midi.track(1).messages(msgCtr).type = 128;
      %midi.track(1).messages(msgCtr).data = [trM(n,3); trM(n,4)];
      %-- note on vel=0:
      midi.track(1).messages(msgCtr).type = 144;
      midi.track(1).messages(msgCtr).data = [trM(n,3); 0];
    end
    msgCtr = msgCtr + 1;
    
    prevtick = cumticks;
  end

  % end of track:
  midi.track(1).messages(msgCtr).deltatime = 0;
  midi.track(1).messages(msgCtr).type = 47;
  midi.track(1).messages(msgCtr).midimeta = 0;
  midi.track(1).messages(msgCtr).data = [];
  midi.track(1).messages(msgCtr).chan = [];
  msgCtr = msgCtr + 1;
  



% return a _column_ vector
% (copied from writemidi.m)
function A=encode_int(val,Nbytes)

A = zeros(Nbytes,1);  %ensure col vector (diff from writemidi.m...)
for i=1:Nbytes
  A(i) = bitand(bitshift(val, -8*(Nbytes-i)), 255);
end

