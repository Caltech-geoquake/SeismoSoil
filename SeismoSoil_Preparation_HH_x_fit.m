function varargout = SeismoSoil_Preparation_HH_x_fit(varargin)
% SEISMOSOIL_PREPARATION_HH_X_FIT MATLAB code for SeismoSoil_Preparation_HH_x_fit.fig
%      SEISMOSOIL_PREPARATION_HH_X_FIT, by itself, creates a new SEISMOSOIL_PREPARATION_HH_X_FIT or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_PREPARATION_HH_X_FIT returns the handle to a new SEISMOSOIL_PREPARATION_HH_X_FIT or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_PREPARATION_HH_X_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_PREPARATION_HH_X_FIT.M with the given input arguments.
%
%      SEISMOSOIL_PREPARATION_HH_X_FIT('Property','Value',...) creates a new SEISMOSOIL_PREPARATION_HH_X_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Preparation_HH_x_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Preparation_HH_x_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Preparation_HH_x_fit

% Last Modified by GUIDE v2.5 08-Sep-2017 00:13:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Preparation_HH_x_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Preparation_HH_x_fit_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Preparation_HH_x_fit is made visible.
function SeismoSoil_Preparation_HH_x_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Preparation_HH_x_fit (see VARARGIN)

% Choose default command line output for SeismoSoil_Preparation_HH_x_fit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Preparation_HH_x_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Preparation_HH_x_fit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_choose_GGmax_data.
function pushbutton1_choose_GGmax_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_choose_GGmax_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select "curve" file...';
[curve_file_name,curve_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');
disp(['User selected ',fullfile(curve_dir_name,curve_file_name)]);
start_dir0 = curve_dir_name;
curve = importdata(fullfile(curve_dir_name,curve_file_name));
handles.metricdata.curve_file_name = curve_file_name;
handles.metricdata.curve_dir_name = curve_dir_name;
handles.metricdata.curve = curve;
handles.metricdata.step1b = 1;

plotCurves(curve);
title(handles.metricdata.curve_file_name,'interpreter','none');

% if handles.metricdata.step1a * handles.metricdata.step1b == 1
if handles.metricdata.step1b == 1  % damping parameter fitting does not require Vs profile
    set(handles.metricdata.handle_Start_damping,'enable','on');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uibuttongroup1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.GGmax_data_source = 'fromVsProfile';
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton1_choose_GGmax_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1_choose_GGmax_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'enable','on');
handle_import_curve = hObject;
handles.metricdata.handle_import_curve = handle_import_curve;
handles.metricdata.step1b = 0;
guidata(hObject,handles);


% --- Executes on button press in pushbutton6_return.
function pushbutton6_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close SeismoSoil_Preparation_HH_x_fit;
SeismoSoil_Preparation;


% --- Executes on button press in checkbox1_show_fig.
function checkbox1_show_fig_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1_show_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1_show_fig

handles.metricdata.show_fig = get(hObject,'value');
guidata(hObject,handles);


% --- Executes on button press in checkbox2_save_fig.
function checkbox2_save_fig_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2_save_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2_save_fig

handles.metricdata.save_fig = get(hObject,'value');
guidata(hObject,handles);


% --- Executes on button press in pushbutton7_damping_para.
function pushbutton7_damping_para_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_damping_para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

save_fig = handles.metricdata.save_fig;
show_fig = handles.metricdata.show_fig;
curve_file_name = handles.metricdata.curve_file_name;
curve_dir_name = handles.metricdata.curve_dir_name;
dir_out = curve_dir_name;

curve_matrix = handles.metricdata.curve;

if strcmpi(curve_file_name(1:5),'curve')
    sitecode = curve_file_name(7:end-4);
else
    sitecode = curve_file_name(1:end-4);
end

if save_fig == 1
    fig_out_dir = uigetdir(start_dir0,'Select folder to save curve-fitting figures...');
    clc;
    [para_xi,fitted_curves] = gaHyX(curve_matrix,show_fig,save_fig,fig_out_dir,sitecode); % "gaHHX" is MKZ + FKZ
else
    clc;
    [para_xi,fitted_curves] = gaHyX(curve_matrix,show_fig,save_fig);  % "gaHyX" = MKZ+muKZ, with "d=1" as the 9th parameter at the end
end

para_xi = [para_xi;ones(1,size(para_xi,2))]; % manually add a row of "1" at the bottom of the matrix

% * * *  Transform para_xi according to the material number  * * * *
try
    vs_profile = handles.metricdata.profile;
    mat = vs_profile(1:end-1,5);
    para_xi_old = para_xi;
    para_xi_new = [];
    for j = 1 : 1 : length(mat)
        para_xi_new = [para_xi_new,para_xi_old(:,mat(j))];
    end
    para_xi = para_xi_new;
catch
    warndlg('You did not choose a Vs profile. Just a reminder...');
end
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

Hy_G_filename = sprintf('HH_X_%s.txt',sitecode);
dlmwrite(fullfile(dir_out,Hy_G_filename),para_xi,'delimiter','\t','precision',6);

fitted_curves_filename = sprintf('Damping_curve_fit_by_HH_%s.txt',sitecode);
dlmwrite(fullfile(dir_out,fitted_curves_filename),fitted_curves,'delimiter','\t','precision',6);


% --- Executes during object creation, after setting all properties.
function pushbutton7_damping_para_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton7_damping_para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'enable','off');
handle_Start_damping = hObject;
handles.metricdata.handle_Start_damping = handle_Start_damping;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uibuttongroup4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.Tmax_data_source = 'fromLadd1991';
handles.metricdata.step2 = 1;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function checkbox1_show_fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox1_show_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.show_fig = 0;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function checkbox2_save_fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox2_save_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.save_fig = 1;
guidata(hObject,handles);


% --- Executes on button press in pushbutton8_close_all.
function pushbutton8_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;
