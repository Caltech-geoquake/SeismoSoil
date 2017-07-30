function varargout = SeismoSoil_Tools_Vs_Profile_Plotter(varargin)
% SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER M-file for SeismoSoil_Tools_Vs_Profile_Plotter.fig
%      SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER, by itself, creates a new SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER returns the handle to a new SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER.M with the given input arguments.
%
%      SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER('Property','Value',...) creates a new SEISMOSOIL_TOOLS_VS_PROFILE_PLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Tools_Vs_Profile_Plotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Tools_Vs_Profile_Plotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Tools_Vs_Profile_Plotter

% Last Modified by GUIDE v2.5 24-Aug-2014 23:45:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Tools_Vs_Profile_Plotter_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Tools_Vs_Profile_Plotter_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Tools_Vs_Profile_Plotter is made visible.
function SeismoSoil_Tools_Vs_Profile_Plotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Tools_Vs_Profile_Plotter (see VARARGIN)

% Choose default command line output for SeismoSoil_Tools_Vs_Profile_Plotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);

% UIWAIT makes SeismoSoil_Tools_Vs_Profile_Plotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Tools_Vs_Profile_Plotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_return_to_tools.
function pushbutton1_return_to_tools_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_return_to_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools_Vs_Profile_Plotter;
SeismoSoil_Tools;


% --- Executes during object creation, after setting all properties.
function pushbutton2_select_profiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2_select_profiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.step4_complete = 0;
guidata(hObject,handles);


% --- Executes on button press in pushbutton2_select_profiles.
function pushbutton2_select_profiles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_select_profiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select Vs profile(s)...';
[motion_file_name,motion_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','on');

if ~isequal(motion_dir_name,0)
    start_dir0 = motion_dir_name;
end

if ischar(motion_file_name) % if motion_file_names is a string
    temp_cell = cell(1,1);
    temp_cell{1} = motion_file_name;
    motion_file_name = cell(1,1);
    motion_file_name{1} = temp_cell{1};
    nr_motion = 1;  % it means that only one motion was selected
else  % otherwise motion_file_names is a cell array
    motion_file_name = motion_file_name.';
    nr_motion = length(motion_file_name);
end

handles.metricdata.motion_file_name = motion_file_name;
handles.metricdata.motion_dir_name = motion_dir_name;
handles.metricdata.nr_motion = nr_motion;

motion = cell(nr_motion,1); % preallocation of cell array
for i = 1 : 1 : nr_motion
    motion{i} = importdata(fullfile(motion_dir_name,motion_file_name{i}));
end
handles.metricdata.motion = motion;
handles.metricdata.step4_complete = 1;

% Initialize baseline correction log
baseline_correction_log = zeros(nr_motion,1);
handles.metricdata.baseline_correction_log = baseline_correction_log;

% Initialize listbox_motions and set selected_motion_indices to 1
handles.metricdata.selected_motion_indices = 1;
temp = handles.metricdata.motion_file_name;
set(handles.listbox_motions,'String',temp,'Value',1);
handles.metricdata.motion_listbox_contents = temp;

guidata(hObject,handles);


% --- Executes on button press in pushbutton3_plot_all.
function pushbutton3_plot_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_plot_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.metricdata.step4_complete == 0
    msgbox('You haven''t selected any motions yet.','Warning');
else

    vs_profile = handles.metricdata.motion;
    nr_profile = handles.metricdata.nr_motion;

    for i = 1 : 1 : nr_profile
        current_profile_filename = handles.metricdata.motion_file_name{i};
        current_profile = vs_profile{i};
        vs30 = calcVs30(current_profile);
        figure_title = sprintf('%s,Vs30=%.1fm/s',current_profile_filename,vs30);
        plotVsProfileForGUI(current_profile,figure_title);
    end
    
end


% --- Executes on button press in pushbutton4_plot_selected.
function pushbutton4_plot_selected_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_plot_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.metricdata.step4_complete == 0
    msgbox('You haven''t selected any motions yet.','Warning');
else
    if handles.metricdata.nr_motion == 1
        selected_motion_indices = 1;
        motion_filenames = handles.metricdata.motion_file_name;
        nr_motion = 1;
    else
        selected_motion_indices = handles.metricdata.selected_motion_indices;
        motion_listbox_contents = handles.metricdata.motion_listbox_contents;
        motion_filenames = motion_listbox_contents(selected_motion_indices);
        nr_motion = length(motion_filenames);
    end

    motion_dir = handles.metricdata.motion_dir_name;
    for i = 1 : 1 : nr_motion
        current_motion_filename = motion_filenames{i};
        current_motion = importdata(fullfile(motion_dir,current_motion_filename));
        plotVsProfileForGUI(current_motion,current_motion_filename);
    end
end


% --- Executes on selection change in listbox_motions.
function listbox_motions_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_motions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_motions

motion_listbox_contents = get(hObject,'String');  
handles.metricdata.motion_listbox_contents = motion_listbox_contents;
selected_motion_indices = get(hObject,'value');
handles.metricdata.selected_motion_indices = selected_motion_indices;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function listbox_motions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5_close_all.
function pushbutton5_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;
