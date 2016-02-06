function varargout = recGUI(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @recGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @recGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function recGUI_OpeningFcn(hObject, eventdata, handles, varargin)
clc;
handles.output = hObject;
set(handles.save,'Enable','off');
set(handles.undo,'visible','off');
guidata(hObject, handles);
global fs;
fs= 44100;
global bits;
bits= 16;
global recObj;
recObj= audiorecorder(fs, bits, 1);
global recordingPath;
recordingPath='rec.wav';
uiwait(handles.figure1);





function varargout = recGUI_OutputFcn(hObject, eventdata, handles) 
global recordingPath;
varargout{1} = recordingPath;
delete(handles.figure1);


function save_Callback(hObject, eventdata, handles)
global fs;
fs1=fs;
global bits
bits1=bits;
global current;
global recordingPath;
plot(current);
set(handles.save,'Enable','off');
%set(handles.text,'string','Your Recording has been saved in Current Folder as "rec.wav".');
[file,path] = uiputfile('*.wav','Save Recording');
fullpath=strcat(path,file);
recordingPath=fullpath;
wavwrite(current, fs1, bits1,recordingPath);
set(handles.text,'string','Your File has been saved.');




function start_Callback(hObject, eventdata, handles)
set(handles.text,'string','');
a=[0];
plot(a);
global fs;
fs1=fs;
global bits
bits1=bits;
global recObj 
% recObj1 = audiorecorder(fs, bits, 1);
record(recObj);



function stop_Callback(hObject, eventdata, handles)
set(handles.text,'string','');
global fs;
fs1=fs;
global bits
bits1=bits;
global recObj; 
% recObj1 = audiorecorder(fs, bits, 1);
stop(recObj);
myRecording = getaudiodata(recObj);
global current;
current=myRecording;
plot(current);
set(handles.save,'Enable','on');
wavwrite(current, fs1, bits1,'rec.wav');


function play_Callback(hObject, eventdata, handles)
set(handles.text,'string','');
global recordingPath;
file=wavread(recordingPath);
wavplay(file,44100);



function clip_Callback(hObject, eventdata, handles)
set(handles.text,'string','');
set(handles.text,'string','Select two Points on the Plot Between which you want to Clip the audio File.');
global fs;
fs1=fs;
global bits
bits1=bits;
global current;
global old;
old=current;
[x1,y1]=ginput(1);
set(handles.text,'string','First Point Recorded. Select Second Point;');
[x2,y2]=ginput(1);
set(handles.text,'string','Second Point Recorded.');
clipped=current(x1:x2);
current=clipped;
set(handles.text,'string','');
plot(current);
set(handles.undo,'visible','on');
set(handles.save,'Enable','on');
wavwrite(current, fs1, bits1,'rec.wav');


function undo_Callback(hObject, eventdata, handles)
set(handles.text,'string','');
global fs;
fs1=fs;
global bits
bits1=bits;
global current;
global old;
current=old;
plot(current);
set(handles.undo,'visible','off');
set(handles.save,'Enable','on');
wavwrite(current, fs1, bits1,'rec.wav');



function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject, 'waitstatus'), 'waiting')
    uiresume(hObject);
    
else
    delete(hObject);
    
end


function pushbutton8_Callback(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);
    
else
    delete(handles.figure1);
    
end
