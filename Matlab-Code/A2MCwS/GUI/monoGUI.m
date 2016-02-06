function varargout = monoGUI(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @monoGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @monoGUI_OutputFcn, ...
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

function monoGUI_OpeningFcn(hObject, eventdata, handles, varargin)

global wav;
global flag;
flag=0;
handles.output = hObject;
set(handles.save,'Enable','off');
set(handles.uipanel2,'visible','off');
set(handles.uipanel1,'visible','off');
set(handles.synth,'Enable','off');
set(handles.listSF,'visible','off');

set(handles.slidHard, 'Min', 0);
set(handles.slidHard, 'Max', 1);
set(handles.slidHard, 'SliderStep', [0.01 0.001]);
set(handles.slidFreq, 'Min', 0);
set(handles.slidFreq, 'Max', 1);
set(handles.slidFreq, 'SliderStep', [0.01 0.001]);

set(handles.thrhval,'string',num2str(get(handles.slidHard,'value')));
set(handles.thrfval,'string',num2str(get(handles.slidFreq,'value')));
guidata(hObject, handles);

function varargout = monoGUI_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;






function GenMIDI_Callback(hObject, eventdata, handles)

global file;
global M;
global Fs;
global flag;
flag=1;


[y,Fs] = wavread(file);

m = 2048; %no of samples in Frame
o = 512; %no of overlapping samples
n = m-o;
w = hanning(m); %hanning window
%w=w/mean(w);
left = y(:,1); %take left channel

time = (1/Fs)*length(left); %Total duration of wave file
t = linspace(0,time,length(left)); %for time axis



thr1=get(handles.slidHard,'value');
left=clip(left,thr1);

[left, k]=padAudio(left,m,n);

thr2=get(handles.slidFreq,'value');
[M,notearray]=detectMono(left,m,n,Fs,w,k,thr2);

M=removeInfFFT(M);

x=str2num(get(handles.edit3,'string'));
M=equalizeNotes(M,n,Fs,x);

M=removeInfFFT(M);
set(handles.uipanel1,'visible','on');
set(handles.uipanel2,'visible','on');
set(handles.listSF,'visible','on');
set(handles.save,'Enable','on');
set(handles.synth,'Enable','on');
global midi_new;
midi_new = matrix2midiHSM(M);
writemidiHSM(midi_new, 'testout.mid');
Notes = midiInfo(midi_new,0);
[PR,t,nn] = piano_roll(Notes,1);
nmat = midi2nmat('testout.mid');
figure (10);
title('Piano Roll');
pianoroll(nmat,'Test','sec','vel');


%%







function synth_Callback(hObject, eventdata, handles)
switch get(get(handles.uipanel1,'SelectedObject'),'string')
    case 'Additive'
        addSynth;
    case 'FM'
        fmSynth;
    case 'Clarinet'
        clarSynth;
end



function save_Callback(hObject, eventdata, handles)

global midi_new;
[file,path] = uiputfile('*.mid','Save Recording');
fullpath=strcat(path,file);

writemidiHSM(midi_new, fullpath);


% --- Executes on button press in add.
function add_Callback(hObject, eventdata, handles)
% hObject    handle to add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add


% --- Executes on button press in clar.
function clar_Callback(hObject, eventdata, handles)
% hObject    handle to clar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of clar


% --- Executes on button press in fm.
function fm_Callback(hObject, eventdata, handles)
% hObject    handle to fm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fm


% --- Executes on selection change in listSF.
function listSF_Callback(hObject, eventdata, handles)
% hObject    handle to listSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listSF contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listSF


% --- Executes during object creation, after setting all properties.
function listSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnPlaySF.
function btnPlaySF_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlaySF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index_selected = get(handles.listSF,'Value');
file_list = get(handles.listSF,'String');
SoundFont = file_list{index_selected};
SoundFont = fullfile(isp_toolboxpath, strcat(SoundFont,'.sf2'));
ispmidi = isp_midiread('testout.mid');
ispmidi.instruments=zeros(1,4);
ispmidi.controller=zeros(1,5);
[wav, fs]=isp_midisynth(ispmidi, 'soundfont', SoundFont);
soundsc(wav,fs);



% --- Executes on button press in btnSaveSF.
function btnSaveSF_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global wav;
index_selected = get(handles.listSF,'Value');
file_list = get(handles.listSF,'String');
SoundFont = file_list{index_selected};
SoundFont = fullfile(isp_toolboxpath, strcat(SoundFont,'.sf2'));
ispmidi = isp_midiread('testout.mid');
ispmidi.instruments=zeros(1,4);
ispmidi.controller=zeros(1,5);
[wav, fs]=isp_midisynth(ispmidi, 'soundfont', SoundFont);
[file,path] = uiputfile('*.wav','Save Recording');
fullpath=strcat(path,file);
wavwrite(wav, fs, 16,fullpath);



function slidHard_Callback(hObject, eventdata, handles)
global flag;
set(handles.thrhval,'string',num2str(get(handles.slidHard,'value')));
if flag==1
    GenMIDI_Callback(hObject, eventdata, handles);
end

function slidFreq_Callback(hObject, eventdata, handles)
global flag;
set(handles.thrfval,'string',num2str(get(handles.slidFreq,'value')));
if flag==1
    GenMIDI_Callback(hObject, eventdata, handles);
end


function slidHard_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function slidFreq_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
