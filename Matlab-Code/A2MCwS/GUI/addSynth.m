function varargout = addSynth(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @addSynth_OpeningFcn, ...
                   'gui_OutputFcn',  @addSynth_OutputFcn, ...
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

function addSynth_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
set(handles.play,'Enable','off');
set(handles.plotAdsr,'Enable','off');
set(handles.saveWav,'Enable','off');
guidata(hObject, handles);


function varargout = addSynth_OutputFcn(hObject, eventdata, handles) 

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
global y;
y=0;
global ADSR;
global Fs;
global At;
global Av;
global Dt;
global Dv;
global St;
% At=At/100;
% Dt=Dt/100;
% St=St/100;
Rt=1-At-Dt-St;
disp(At);
disp(Av);
disp(Dt);
disp(Dv);
disp(St);
if (get(handles.togglebutton1,'value'))==0
    b=str2num(get(handles.edit1,'string'));
    c=str2num(get(handles.edit3,'string'));
    wavef1=get(get(handles.uipanel4,'SelectedObject'),'string');
    wavef2=get(get(handles.uipanel5,'SelectedObject'),'string');
    wavef3=get(get(handles.uipanel6,'SelectedObject'),'string');
    x1=str2func(wavef1);
    x2=str2func(wavef2);
    x3=str2func(wavef3);
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
         c=a*c;
         y = [y (a*x1(f*2*pi*t)+b*x2(f*pi*t)+c*x3(f*4*pi*t)).*ADSR(1,:)];
    end
         else
             for i=1:length(M(:,1))
                t=M(i,5):1/Fs:M(i,6);
                f= (2^((M(i,3)-69)/12))*(440);
                Ae=linspace(0,Av,At*length(t));
                De=linspace(Av,Dv,Dt*length(t));
                Se=linspace(Dv,Dv,St*length(t));
                Re=linspace(Dv,0,Rt*length(t));
                ADSR=[Ae De Se Re];
                ADSR=padarray(ADSR,[0 (length(t)-length(ADSR))],0,'post');
                w=2*pi*f*t;
                a=M(i,4)/127;
                y = [y (a*(sin(w)+0.75*sin(3*w)+0.5*sin(5*w)+0.14*sin(7*w)+0.5*sin(9*w)+0.14*sin(11*w)+0.17*sin(13*w))).*ADSR(1,:)];
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



function togglebutton1_Callback(hObject, eventdata, handles)
set(handles.play,'Enable','off');
set(handles.plotAdsr,'Enable','off');
set(handles.saveWav,'Enable','off');
if (get(handles.togglebutton1,'value'))==1
    set(handles.togglebutton1,'string','Click for Additive Synthesis');
    set(handles.uipanel7,'visible','off');
else 
    set(handles.uipanel7,'visible','on');
end
