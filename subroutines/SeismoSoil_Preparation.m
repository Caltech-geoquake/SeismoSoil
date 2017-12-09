function varargout = SeismoSoil_Preparation(varargin)
% SEISMOSOIL_PREPARATION MATLAB code for SeismoSoil_Preparation.fig
%      SEISMOSOIL_PREPARATION, by itself, creates a new SEISMOSOIL_PREPARATION or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_PREPARATION returns the handle to a new SEISMOSOIL_PREPARATION or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_PREPARATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_PREPARATION.M with the given input arguments.
%
%      SEISMOSOIL_PREPARATION('Property','Value',...) creates a new SEISMOSOIL_PREPARATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Preparation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Preparation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Preparation

% Last Modified by GUIDE v2.5 08-Sep-2017 00:16:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Preparation_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Preparation_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Preparation is made visible.
function SeismoSoil_Preparation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Preparation (see VARARGIN)

% Choose default command line output for SeismoSoil_Preparation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Preparation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Preparation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_Vs_profile.
function pushbutton1_Vs_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_Vs_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Preparation;
SeismoSoil_Preparation_Vs_Profile;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Preparation;
SeismoSoil_Preparation_Modulus_Damping;


% --- Executes on button press in pushbutton5_curve_fitting.
function pushbutton5_curve_fitting_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_curve_fitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Preparation;
SeismoSoil_Preparation_H2_H4_Curve_Fitting;


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Preparation;
SeismoSoil_Preparation_HH_G_fit;


% --- Executes on button press in pushbutton7_HH_x_fit.
function pushbutton7_HH_x_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_HH_x_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Preparation;
SeismoSoil_Preparation_HH_x_fit;


% --- Executes on button press in pushbutton4_return.
function pushbutton4_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close SeismoSoil_Preparation;
SeismoSoil;
