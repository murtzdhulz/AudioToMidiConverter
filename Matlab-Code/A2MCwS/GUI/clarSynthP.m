function varargout = clarSynthP(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @clarSynthP_OpeningFcn, ...
                   'gui_OutputFcn',  @clarSynthP_OutputFcn, ...
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

function clarSynthP_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
set(handles.play,'Enable','off');
set(handles.plotAdsr,'Enable','off');
set(handles.saveWav,'Enable','off');

guidata(hObject, handles);

function varargout = clarSynthP_OutputFcn(hObject, eventdata, handles) 

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



function synth_Callback(hObject, eventdata, handles)
global M;
global Fs;
global At;
global Av;
global Dt;
global Dv;
global St;
Rt=1-At-Dt-St;
global y;
global ADSR;
y=0;
tt = M(1,5):1/Fs:M(length(M(:,6)),6);
y = zeros(1,length(tt));
for i=1:length(M(:,1))
    t=M(i,5):1/Fs:M(i,6);
    ind=round(M(i,5)*Fs);
    f= (2^((M(i,3)-69)/12))*(440);
    Ae=linspace(0,Av,At*length(t))
    De=linspace(Av,Dv,Dt*length(t));
    Se=linspace(Dv,Dv,St*length(t));
    Re=linspace(Dv,0,Rt*length(t));
    ADSR=[Ae De Se Re];
    ADSR=padarray(ADSR,[0 (length(t)-length(ADSR))],0,'post');
    w=2*pi*f*t;
    a=M(i,4)/127;
    z=(a*(sin(w)+0.75*sin(3*w)+0.5*sin(5*w)+0.14*sin(7*w)+0.5*sin(9*w)+0.14*sin(11*w)+0.17*sin(13*w))).*ADSR(1,:);
    for j=1:length(t)
        y(j+ind) = y(j+ind)+z(j);
    end
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
