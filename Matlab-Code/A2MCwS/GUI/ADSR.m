function varargout = ADSR(varargin)
% ADSR MATLAB code for ADSR.fig
%      ADSR, by itself, creates a new ADSR or raises the existing
%      singleton*.
%
%      H = ADSR returns the handle to a new ADSR or the handle to
%      the existing singleton*.
%
%      ADSR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADSR.M with the given input arguments.
%
%      ADSR('Property','Value',...) creates a new ADSR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ADSR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ADSR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ADSR

% Last Modified by GUIDE v2.5 15-Apr-2014 22:10:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ADSR_OpeningFcn, ...
                   'gui_OutputFcn',  @ADSR_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before ADSR is made visible.
function ADSR_OpeningFcn(hObject, eventdata, handles, varargin)
global At;
global Av;
global Dt;
global Dv;
global St;
handles.output = hObject;
plot(1);
axis([0 100 0 1]);
guidata(hObject, handles);
set(handles.text1,'string','Click to select attack time and attack value');
[At,Av]=ginput(1);
if Av>1
    Av=1;
end
x=[0,At];
y=[0,Av];
plot(x,y);
axis([0 100 0 1]);

set(handles.text1,'string','Click to select decay time and decay value');
[Dt,Dv]=ginput(1);
Dt=Dt-At;
x=[0,At,Dt+At];
y=[0,Av,Dv];
plot(x,y);
axis([0 100 0 1]);

set(handles.text1,'string','Click to select sustain time');
[St,a]=ginput(1);
St=St-Dt-At;
x=[0,At,Dt+At,St+Dt+At];
y=[0,Av,Dv,Dv];
plot(x,y);
axis([0 100 0 1]);

set(handles.text1,'string','');
x=[0,At,Dt+At,St+Dt+At,100];
y=[0,Av,Dv,Dv,0];
plot(x,y);
axis([0 100 0 1]);

At=At/100;
Dt=Dt/100;
St=St/100;

% UIWAIT makes ADSR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ADSR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
