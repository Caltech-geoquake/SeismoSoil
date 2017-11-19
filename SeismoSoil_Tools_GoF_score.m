function varargout = SeismoSoil_Tools_GoF_score(varargin)
% SEISMOSOIL_TOOLS_GOF_SCORE MATLAB code for SeismoSoil_Tools_GoF_score.fig
%      SEISMOSOIL_TOOLS_GOF_SCORE, by itself, creates a new SEISMOSOIL_TOOLS_GOF_SCORE or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_TOOLS_GOF_SCORE returns the handle to a new SEISMOSOIL_TOOLS_GOF_SCORE or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_TOOLS_GOF_SCORE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_TOOLS_GOF_SCORE.M with the given input arguments.
%
%      SEISMOSOIL_TOOLS_GOF_SCORE('Property','Value',...) creates a new SEISMOSOIL_TOOLS_GOF_SCORE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Tools_GoF_score_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Tools_GoF_score_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Tools_GoF_score

% Last Modified by GUIDE v2.5 04-Sep-2017 02:12:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Tools_GoF_score_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Tools_GoF_score_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Tools_GoF_score is made visible.
function SeismoSoil_Tools_GoF_score_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Tools_GoF_score (see VARARGIN)

% Choose default command line output for SeismoSoil_Tools_GoF_score
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Tools_GoF_score wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%cwd = pwd();  % get current (absolute) working directory
%addpath(fullfile(cwd,'lib'));  % let MATLAB search in ./lib for all necessary subroutines

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Tools_GoF_score_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_select_meas.
function pushbutton_select_meas_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_meas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select measurement data file...';
[meas_file_name,meas_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');

start_dir0 = meas_dir_name;

handles.metricdata.meas_file_name = meas_file_name;
handles.metricdata.meas_dir_name = meas_dir_name;

handles.metricdata.step1_complete = 1;

set(handles.edit1_path_meas,'string',sprintf('%s...%s%s',meas_dir_name(1:21),meas_dir_name(end-4:end),meas_file_name));

guidata(hObject,handles);


% --- Executes on button press in pushbutton_select_simu.
function pushbutton_select_simu_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_simu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select simulation data file...';
[simu_file_name,simu_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');

start_dir0 = simu_dir_name;

handles.metricdata.simu_file_name = simu_file_name;
handles.metricdata.simu_dir_name = simu_dir_name;

handles.metricdata.step2_complete = 1;

set(handles.edit2_path_simu,'string',sprintf('%s...%s%s',simu_dir_name(1:21),simu_dir_name(end-4:end),simu_file_name));

guidata(hObject,handles);



function edit1_path_meas_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_path_meas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_path_meas as text
%        str2double(get(hObject,'String')) returns contents of edit1_path_meas as a double

meas_full_filename = get(hObject,'String');
[meas_dir_name,meas_file_name,meas_ext] = fileparts(meas_full_filename);
handles.metricdata.meas_dir_name = meas_dir_name;
handles.metricdata.meas_file_name = sprintf('%s%s',meas_file_name,meas_ext);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit1_path_meas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_path_meas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_path_simu_Callback(hObject, eventdata, handles)
% hObject    handle to edit2_path_simu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2_path_simu as text
%        str2double(get(hObject,'String')) returns contents of edit2_path_simu as a double

simu_full_filename = get(hObject,'String');
[simu_dir_name,simu_file_name,simu_ext] = fileparts(simu_full_filename);
handles.metricdata.simu_dir_name = simu_dir_name;
handles.metricdata.simu_file_name = sprintf('%s%s',simu_file_name,simu_ext);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit2_path_simu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2_path_simu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_fmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_fmin as text
%        str2double(get(hObject,'String')) returns contents of edit3_fmin as a double

handles.metricdata.fmin = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit3_fmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_fmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit4_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4_fmax as text
%        str2double(get(hObject,'String')) returns contents of edit4_fmax as a double

handles.metricdata.fmax = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit4_fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.meas_unit = 'm/s/s';
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1a'
        meas_unit = 'm/s/s';
    case 'radiobutton1b'
        meas_unit = 'gal';
    case 'radiobutton1c'
        meas_unit = 'g';
end
handles.metricdata.meas_unit = meas_unit;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.simu_unit = 'm/s/s';
guidata(hObject,handles);



% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton2a'
        simu_unit = 'm/s/s';
    case 'radiobutton2b'
        simu_unit = 'gal';
    case 'radiobutton2c'
        simu_unit = 'g';
end
handles.metricdata.simu_unit = simu_unit;
guidata(hObject,handles);


% --- Executes on button press in pushbutton3_calc_score.
function pushbutton3_calc_score_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_calc_score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

meas_dir = handles.metricdata.meas_dir_name;
meas_file = handles.metricdata.meas_file_name;

simu_dir = handles.metricdata.simu_dir_name;
simu_file = handles.metricdata.simu_file_name;

meas = importdata(fullfile(meas_dir,meas_file));
simu = importdata(fullfile(simu_dir,simu_file));

switch handles.metricdata.meas_unit
    case 'gal'
        meas(:,2) = meas(:,2) / 100;
    case 'g'
        meas(:,2) = meas(:,2) * 9.81;
end
   
switch handles.metricdata.simu_unit
    case 'gal'
        simu(:,2) = simu(:,2) / 100;
    case 'g'
        simu(:,2) = simu(:,2) * 9.81;
end

dt_meas = meas(2,1) - meas(1,1);
dt_simu = simu(2,1) - simu(1,1);

if dt_meas > dt_simu
    meas_ = meas;
    simu_ = [meas(:,1), interp1(simu(:,1),simu(:,2),meas(:,1))];  % downsample
elseif dt_meas < dt_simu
    simu_ = simu;
    meas_ = [simu(:,1), interp1(meas(:,1),meas(:,2),simu(:,1))];  % downsample
else
    meas_ = meas;
    simu_ = simu;
end

fmin = handles.metricdata.fmin;
fmax = handles.metricdata.fmax;

if fmin <= 0
    fprintf('****  Warning: f_min value is invalid; using 0.1 Hz instead.   ***\n');
    fmin = 0.1;
end

fprintf('Calculating goodness-of-fit scores...  \n\n');
[scores,avg_score] = gofScores(meas_,simu_,fmin,fmax,[1 1 1 1],[1 1],4,'n');
fprintf('Individual scores:\n')
fprintf('  S1 = %5.1f, S2 = %5.1f\n',scores(1),scores(2));
fprintf('  S3 = %5.1f, S4 = %5.1f\n',scores(3),scores(4));
fprintf('  S5 = %5.1f, S6 = %5.1f\n',scores(5),scores(6));
fprintf('  S7 = %5.1f, S8 = %5.1f\n',scores(7),scores(8));
fprintf('  S9 = %5.1f\n',scores(9));
fprintf('Average score = %.1f\n',avg_score);


% --- Executes on button press in pushbutton4_close_all.
function pushbutton4_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;
close SeismoSoil_Tools_GoF_score;
SeismoSoil_Tools;
