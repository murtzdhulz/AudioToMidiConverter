function varargout = fileipGUI(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fileipGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @fileipGUI_OutputFcn, ...
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

function fileipGUI_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
clc;
guidata(hObject, handles);
set(handles.path,'string','Selected File : ');
set(handles.undo,'visible','off');
set(handles.save,'Enable','off');
global fs;
fs= 44100;
global bits;
bits= 16;
global inputPath;
inputPath='rec.wav';
uiwait(handles.figure1);


function varargout = fileipGUI_OutputFcn(hObject, eventdata, handles) 
global inputPath;
varargout{1} = inputPath;
delete(handles.figure1);



function browse_Callback(hObject, eventdata, handles)
[file1 path]=uigetfile({'*.wav'},'Select a WAVE File');
global inputPath;
fullpath=strcat(path,file1);
inputPath=fullpath;
global file;
file=wavread(fullpath);
plot(file);
set(handles.text,'string',fullpath);
global current;
current=file;


function play_Callback(hObject, eventdata, handles)
global current;
wavplay(current,44100);



function clip_Callback(hObject, eventdata, handles)
set(handles.textb,'string','');
set(handles.textb,'string','Select two Points Between on the Plot which you want to Clip the audio File.');
global current;
global old;
old=current;
[x1,y1]=ginput(1);
set(handles.textb,'string','First Point Recorded. Select Second Point;');
[x2,y2]=ginput(1);
set(handles.textb,'string','Second Point Recorded.');
clipped=current(x1:x2);
current=clipped;
set(handles.text,'string','');
plot(current);
set(handles.undo,'visible','on');
set(handles.save,'Enable','on');

function undo_Callback(hObject, eventdata, handles)
set(handles.textb,'string','');
global current;
global old;
current=old;
plot(current);
set(handles.undo,'visible','off');
set(handles.save,'Enable','on');



function save_Callback(hObject, eventdata, handles)
set(handles.textb,'string','');
global fs;
fs1=fs;
global bits
bits1=bits;
global current;
global inputPath;
[file,path] = uiputfile('*.wav','Save Recording');
fullpath=strcat(path,file);
inputPath=fullpath;
wavwrite(current, fs1, bits1,inputPath);
set(handles.textb,'string','Your File has been saved.');



function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject, 'waitstatus'), 'waiting')
    uiresume(hObject);
else
    delete(hObject);
end


function pushbutton15_Callback(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);
else
    delete(handles.figure1);
end
