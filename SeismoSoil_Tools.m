function varargout = SeismoSoil_Tools(varargin)
% SEISMOSOIL_TOOLS M-file for SeismoSoil_Tools.fig
%      SEISMOSOIL_TOOLS, by itself, creates a new SEISMOSOIL_TOOLS or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_TOOLS returns the handle to a new SEISMOSOIL_TOOLS or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_TOOLS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_TOOLS.M with the given input arguments.
%
%      SEISMOSOIL_TOOLS('Property','Value',...) creates a new SEISMOSOIL_TOOLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Tools_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Tools_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Tools

% Last Modified by GUIDE v2.5 21-May-2015 15:18:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Tools_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Tools_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Tools is made visible.
function SeismoSoil_Tools_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Tools (see VARARGIN)

% Choose default command line output for SeismoSoil_Tools
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Tools wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Tools_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in pushbutton1_return_to_main.
function pushbutton1_return_to_main_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_return_to_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil;


% --- Executes on button press in pushbutton2_motion_plotter.
function pushbutton2_motion_plotter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_motion_plotter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_Motion_Plotter;


% --- Executes on button press in pushbutton3_fourier_transform.
function pushbutton3_fourier_transform_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_fourier_transform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_Fourier_Transform;


% --- Executes on button press in pushbutton4_response_spectra.
function pushbutton4_response_spectra_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_response_spectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_Response_Spectra;


% --- Executes on button press in pushbutton5_Vs_profile_plotter.
function pushbutton5_Vs_profile_plotter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_Vs_profile_plotter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_Vs_Profile_Plotter;


% --- Executes on button press in pushbutton6_baseline_correction.
function pushbutton6_baseline_correction_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6_baseline_correction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_Baseline_Correction;


% --- Executes on button press in pushbutton7_filter.
function pushbutton7_filter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_Signal_Filters;


% --- Executes on button press in pushbutton8_deconvolution.
function pushbutton8_deconvolution_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8_deconvolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_Deconvolution;


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_Motion_Unit_Converter;


% --- Executes on button press in pushbutton10_PEEG_to_2col.
function pushbutton10_PEEG_to_2col_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10_PEEG_to_2col (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_PEER_To_2COL;


% --- Executes on button press in pushbutton11_SMC_to_2col.
function pushbutton11_SMC_to_2col_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11_SMC_to_2col (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools;
SeismoSoil_Tools_SMC_To_2COL;
