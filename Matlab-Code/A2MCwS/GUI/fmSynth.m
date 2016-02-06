function varargout = fmSynth(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fmSynth_OpeningFcn, ...
                   'gui_OutputFcn',  @fmSynth_OutputFcn, ...
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

function fmSynth_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
set(handles.play,'Enable','off');
set(handles.plotAdsr,'Enable','off');
set(handles.saveWav,'Enable','off');

guidata(hObject, handles);

function varargout = fmSynth_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)

function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function synth_Callback(hObject, eventdata, handles)
global M;
global Fs;
b=str2num(get(handles.edit1,'string'));
fc=str2num(get(handles.edit2,'string'));
wavef1=get(get(handles.uipanel2,'SelectedObject'),'string');
wavef2=get(get(handles.uipanel1,'SelectedObject'),'string');
global At;
global Av;
global Dt;
global Dv;
global St;
Rt=1-At-Dt-St;
global y;
global ADSR;
y=0;
    x1=str2func(wavef1);
    x2=str2func(wavef2);
    for i=1:length(M(:,1))
        t=M(i,5):1/Fs:M(i,6);
        f= (2^((M(i,3)-69)/12))*(440);
        Ae=linspace(0,Av,At*length(t));
        De=linspace(Av,Dv,Dt*length(t));
        Se=linspace(Dv,Dv,St*length(t));
        Re=linspace(Dv,0,Rt*length(t));
        ADSR=[Ae De Se Re];
        ADSR=padarray(ADSR,[0 (length(t)-length(ADSR))],0,'post');
        a=M(i,4)/127;
        b=a*b;
        y = [y (b*x2(2*fc*pi*t+a*x1(f*2*pi*t))).*ADSR(1,:)];
    end
set(handles.play,'Enable','on');
set(handles.plotAdsr,'Enable','on');
set(handles.saveWav,'Enable','on');


function plotAdsr_Callback(hObject, eventdata, handles)
global ADSR;
figure;
plot(ADSR);

function play_Callback(hObject, eventdata, handles)
global y;
global Fs;
soundsc(y,Fs);

function saveWav_Callback(hObject, eventdata, handles)
global y;
global Fs;
[file,path] = uiputfile('*.wav','Save Recording');
fullpath=strcat(path,file);
wavwrite(y, Fs, 16,fullpath);
