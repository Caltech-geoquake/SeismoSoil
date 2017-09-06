function varargout = SeismoSoil(varargin)
% SeismoSoil M-file for SeismoSoil.fig
%      SeismoSoil, by itself, creates a new SeismoSoil or raises the existing
%      singleton*.
%
%      H = SeismoSoil returns the handle to a new SeismoSoil or the handle to
%      the existing singleton*.
%
%      SeismoSoil('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SeismoSoil.M with the given input arguments.
%
%      SeismoSoil('Property','Value',...) creates a new SeismoSoil or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil

% Last Modified by GUIDE v2.5 20-May-2015 13:26:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_OutputFcn, ...
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


% --- Executes just before SeismoSoil is made visible.
function SeismoSoil_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil (see VARARGIN)

% Choose default command line output for SeismoSoil
handles.output = hObject;

% Global "uigetfile start dir"
% % handles.metricdata.uigetfile_start_dir = pwd;
global start_dir0;
if isempty(start_dir0)
    start_dir0 = pwd();
end

[SS_dir,~,~] = fileparts(mfilename('fullpath'));  % (absolute) directory where SeismoSoil.m is stored
addpath(fullfile(SS_dir,'lib'));  % let MATLAB search in ./lib for all necessary subroutines

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton2_equivalent_linear.
function pushbutton2_equivalent_linear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_equivalent_linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Equivalent_Linear_Analysis;



% --- Executes on button press in pushbutton3_nonlinear_H2.
function pushbutton3_nonlinear_H2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_nonlinear_H2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Nonlinear_Analysis_H2;


% --- Executes on button press in pushbutton1_input_preparation.
function pushbutton1_input_preparation_Callback(hObject,eventdata,handles)
% hObject    handle to pushbutton1_input_preparation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Input_Files_Preperation;

% % % --------------------------------------------------------------------
% % function Menu1_help_Callback(hObject, eventdata, handles)
% % % hObject    handle to Menu1_help (see GCBO)
% % % eventdata  reserved - to be defined in a future version of MATLAB
% % % handles    structure with handles and user data (see GUIDATA)
% % 
% % 
% % 
% % % --------------------------------------------------------------------
% % function Menu1_1_about_Callback(hObject, eventdata, handles)
% % % hObject    handle to Menu1_1_about (see GCBO)
% % % eventdata  reserved - to be defined in a future version of MATLAB
% % % handles    structure with handles and user data (see GUIDATA)
% % 
% % msgbox({'v1.1.4, 9/21/2013','(First draft: 6/30/2013)',' ','Authors: Dominic Assimaki, Wei Li & Jian Shi',...
% %     'GUI Design: Jian Shi',' ','Contact: Dr. Assimaki - dominic.assimaki@ce.gatech.edu',...
% %     '              Jian Shi - shijian@gatech.edu'},'About');
% % 


% --- Executes on button press in pushbutton4_tools.
function pushbutton4_tools_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Tools;


% --- Executes on button press in pushbutton5_linear_analysis.
function pushbutton5_linear_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_linear_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Linear_Analysis;


% --- Executes on button press in pushbutton6_nonlinear_H4.
function pushbutton6_nonlinear_H4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6_nonlinear_H4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Nonlinear_Analysis_H4;


% --- Executes on button press in pushbutton12_nonlinear_Hy.
function pushbutton12_nonlinear_Hy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12_nonlinear_Hy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil
SeismoSoil_Nonlinear_Analysis_Hy;


% --- Executes on button press in pushbutton14_nonlinear_HH.
function pushbutton14_nonlinear_HH_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14_nonlinear_HH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil
SeismoSoil_Nonlinear_Analysis_HH;


% --- Executes on button press in pushbutton7_fewer_clicks.
function pushbutton7_fewer_clicks_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_fewer_clicks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close SeismoSoil;
SeismoSoil_Fewer_Clicks_Analysis;


% --- Executes on button press in pushbutton8_Linear_Time_Domain.
function pushbutton8_Linear_Time_Domain_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8_Linear_Time_Domain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Linear_Analysis_Time_Domain;


% --- Executes on button press in pushbutton13_nonlinear_EPP.
function pushbutton13_nonlinear_EPP_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13_nonlinear_EPP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Nonlinear_Analysis_EPP;


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil;
SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis;


% --- Executes on button press in pushbutton_author.
function pushbutton_author_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_author (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hm = msgbox({'Current version: v1.2.8.5, 9/5/2017',...
    'First draft: 6/30/2013',...
    ' ','Authors: Domniki Asimaki, Jian Shi, Wei Li',...
    'GUI Design: Jian Shi',' ','Contact: D.A.- domniki@caltech.edu',...
    '             J.Sh.- jian.shi@caltech.edu'},'About');
ahm = get(hm,'CurrentAxes');
chm = get(ahm,'Children');
set(chm,'FontSize',10);
set(ahm,'FontSize',10);
set(chm,'FontName','Helvetica');
pos = get(hm,'position');
pos = pos + [0,0,42,10];
set(hm,'position',pos);

% Change log:
% 03/18/2014 - Added PI (plasticity index) in GTSRA_Dynamic_Soil_Properties
% 06/07/2014 - Various updats of signal filtering, Fourier spectra, baseline
%              correction, etc.
% 06/08/2014 - Bug fix in fourierTransform


