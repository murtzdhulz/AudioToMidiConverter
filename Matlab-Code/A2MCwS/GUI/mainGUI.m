function varargout = mainGUI(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mainGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mainGUI_OutputFcn, ...
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

function mainGUI_OpeningFcn(hObject, eventdata, handles, varargin)
global file;

handles.output = hObject;


guidata(hObject, handles);



function varargout = mainGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function select_Callback(hObject, eventdata, handles)
global file;
file=fileipGUI;
set(handles.path,'string',file);

function record_Callback(hObject, eventdata, handles)
global file;
file=recGUI;
set(handles.path,'string',file);



function mono_Callback(hObject, eventdata, handles)
global file;
h=handles.figure1;
setappdata(h,'FilePath',file);
monoGUI;




function poly_Callback(hObject, eventdata, handles)
global file;
h=handles.figure1;
setappdata(h,'FilePath',file);
polyGUI;
    



function adsr_Callback(hObject, eventdata, handles)
ADSR;
