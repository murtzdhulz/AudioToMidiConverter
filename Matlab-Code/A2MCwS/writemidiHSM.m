function rawbytes=writemidi(midi,filename,do_run_mode)

do_run_mode = 0;

databytes_track{1} = [];
  
  for j=1:length(midi.track(1).messages)

    msg = midi.track(1).messages(j);

    msg_bytes = encode_var_length(msg.deltatime);

    if (msg.midimeta==1)

      % check for doing running mode
      run_mode = 0;
      run_mode = msg.used_running_mode;
      
      % should check that prev msg has same type to allow run
      % mode...
      
      
      %      if (j>1 && do_run_mode && msg.type == midi.track(1).messages(j-1).type)
%	run_mode = 1;
%      end


msg_bytes = [msg_bytes; encode_midi_msg(msg, run_mode)];
    
    
    else
      
      msg_bytes = [msg_bytes; encode_meta_msg(msg)];
      
    end



    databytes_track{1} = [databytes_track{1}; msg_bytes];
    
  end
  
  % HEADER
% double('MThd') = [77 84 104 100]
rawbytes = [77 84 104 100 ...
	    0 0 0 6 ...
	    encode_int(midi.format,2) ...
	    encode_int(1,2) ...
	    encode_int(midi.ticks_per_quarter_note,2) ...
	   ]';

% TRACK_CHUCKS

  a = length(databytes_track{1});
  % double('MTrk') = [77 84 114 107]
  tmp = [77 84 114 107 ...
	 encode_int(length(databytes_track{1}),4) ...
	 databytes_track{1}']';
  rawbytes(end+1:end+length(tmp)) = tmp;
  
  % write to file
fid = fopen(filename,'w');
%fwrite(fid,rawbytes,'char');
fwrite(fid,rawbytes,'uint8');
fclose(fid);

% return a _column_ vector
function A=encode_int(val,Nbytes)

for i=1:Nbytes
  A(i) = bitand(bitshift(val, -8*(Nbytes-i)), 255);
end


function bytes=encode_var_length(val)
if val<128 % covers 99% cases!
    bytes = val;
    return
end
binStr = dec2base(round(val),2);
Nbytes = ceil(length(binStr)/7);
binStr = ['00000000' binStr];
bytes = [];
for i=1:Nbytes
  if (i==1)
    lastbit = '0';
  else
    lastbit = '1';
  end
  B = bin2dec([lastbit binStr(end-i*7+1:end-(i-1)*7)]);
  bytes = [B; bytes];
end


function bytes=encode_midi_msg(msg, run_mode)

bytes = [];

if (run_mode ~= 1)
  bytes = msg.type;
  % channel:
  bytes = bytes + msg.chan;  % lower nibble should be chan
end

bytes = [bytes; msg.data];

function bytes=encode_meta_msg(msg)

bytes = 255;
bytes = [bytes; msg.type];
bytes = [bytes; encode_var_length(length(msg.data))];
bytes = [bytes; msg.data];


  
  