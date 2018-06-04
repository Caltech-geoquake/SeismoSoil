function varargout = SeismoSoil_Tools_Response_Spectra(varargin)
% SEISMOSOIL_TOOLS_RESPONSE_SPECTRA M-file for SeismoSoil_Tools_Response_Spectra.fig
%      SEISMOSOIL_TOOLS_RESPONSE_SPECTRA, by itself, creates a new SEISMOSOIL_TOOLS_RESPONSE_SPECTRA or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_TOOLS_RESPONSE_SPECTRA returns the handle to a new SEISMOSOIL_TOOLS_RESPONSE_SPECTRA or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_TOOLS_RESPONSE_SPECTRA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_TOOLS_RESPONSE_SPECTRA.M with the given input arguments.
%
%      SEISMOSOIL_TOOLS_RESPONSE_SPECTRA('Property','Value',...) creates a new SEISMOSOIL_TOOLS_RESPONSE_SPECTRA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Tools_Response_Spectra_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Tools_Response_Spectra_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Tools_Response_Spectra

% Last Modified by GUIDE v2.5 24-Aug-2014 23:45:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Tools_Response_Spectra_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Tools_Response_Spectra_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Tools_Response_Spectra is made visible.
function SeismoSoil_Tools_Response_Spectra_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Tools_Response_Spectra (see VARARGIN)

% Choose default command line output for SeismoSoil_Tools_Response_Spectra
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);

% UIWAIT makes SeismoSoil_Tools_Response_Spectra wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Tools_Response_Spectra_OutputFcn(hObject, eventdata, handles) 
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

close SeismoSoil_Tools_Response_Spectra;
SeismoSoil_Tools;

% --- Executes during object creation, after setting all properties.
function pushbutton2_select_motions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2_select_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.step4_complete = 0;
guidata(hObject,handles);


% --- Executes on button press in pushbutton2_select_motions.
function pushbutton2_select_motions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_select_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select input motion...';
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
    checkInputs(motion{i},'motion',sprintf('Motion #%d',i));
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


% % % 
% % % % --- Executes on button press in pushbutton2.
% % % function pushbutton2_Callback(hObject, eventdata, handles)
% % % % hObject    handle to pushbutton2 (see GCBO)
% % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % handles    structure with handles and user data (see GUIDATA)
% % % 
% % % 
% % % % --- Executes on button press in pushbutton3.
% % % function pushbutton3_Callback(hObject, eventdata, handles)
% % % % hObject    handle to pushbutton3 (see GCBO)
% % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % handles    structure with handles and user data (see GUIDATA)


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

% % % 
% % % % --- Executes on button press in pushbutton4.
% % % function pushbutton4_Callback(hObject, eventdata, handles)
% % % % hObject    handle to pushbutton4 (see GCBO)
% % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function uipanel1_x_axis_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1_x_axis_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% % accel_unit_division_factor = 1;  % default value of unit selection
% % input_accel_unit_ = 'm/s^2'; % default value of unit selection
% % handles.metricdata.accel_unit_division_factor = accel_unit_division_factor;
% % handles.metricdata.input_accel_unit = input_accel_unit_;

x_axis_scale = 'log';
handles.metricdata.x_axis_scale = x_axis_scale;
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel1_x_axis_scale.
function uipanel1_x_axis_scale_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1_x_axis_scale 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1a'
        x_axis_scale = 'linear';
    case 'radiobutton1b'
        x_axis_scale = 'log';
end
handles.metricdata.x_axis_scale = x_axis_scale;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function uipanel8_unit_of_accel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel8_unit_of_accel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

accel_unit = 'm/s^2';
handles.metricdata.accel_unit = accel_unit;
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel8_unit_of_accel.
function uipanel8_unit_of_accel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel8_unit_of_accel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton2a'
        accel_unit = 'm/s^2';
    case 'radiobutton2b'
        accel_unit = 'gal';
    case 'radiobutton2c'
        accel_unit = 'g';
end
handles.metricdata.accel_unit = accel_unit;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit1_damping_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_damping_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

damping_ratio = 0.05;
handles.metricdata.damping_ratio = damping_ratio;
guidata(hObject,handles);



function edit1_damping_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_damping_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_damping_ratio as text
%        str2double(get(hObject,'String')) returns contents of edit1_damping_ratio as a double

damping_ratio = str2double(get(hObject,'String'))/100; % convert from percent to 1
handles.metricdata.damping_ratio = damping_ratio;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit2_array_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2_array_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

array_density = 100;
handles.metricdata.array_density = array_density;
guidata(hObject,handles);


function edit2_array_density_Callback(hObject, eventdata, handles)
% hObject    handle to edit2_array_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2_array_density as text
%        str2double(get(hObject,'String')) returns contents of edit2_array_density as a double

array_density = str2double(get(hObject,'String'));
handles.metricdata.array_density = array_density;
guidata(hObject,handles);


% --- Executes on button press in pushbutton7_whats_this.
function pushbutton7_whats_this_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_whats_this (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox({'This fields defines how many uniformly-distributed',...
        'points there will be in an order of magnitude',...
        '(i.e., from 1 to 10, from 100 to 1000, etc.).'},...
        'Help');

% --- Executes during object creation, after setting all properties.
function edit3_maximum_period_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_maximum_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

max_period = 10;
handles.metricdata.max_period = max_period;
guidata(hObject,handles);

    
function edit3_maximum_period_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_maximum_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_maximum_period as text
%        str2double(get(hObject,'String')) returns contents of edit3_maximum_period as a double

max_period = str2double(get(hObject,'String'));
handles.metricdata.max_period = max_period;
guidata(hObject,handles);




% % % % --- Executes during object creation, after setting all properties.
% % % function uipanel2_double_sided_CreateFcn(hObject, eventdata, handles)
% % % % hObject    handle to uipanel2_double_sided (see GCBO)
% % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % handles    empty - handles not created until after all CreateFcns called
% % % 
% % % side_option = 'single';
% % % handles.metricdata.side_option = side_option;
% % % guidata(hObject,handles);
% % % 
% % % 
% % % 
% % % % --- Executes when selected object is changed in uipanel2_double_sided.
% % % function uipanel2_double_sided_SelectionChangeFcn(hObject, eventdata, handles)
% % % % hObject    handle to the selected object in uipanel2_double_sided 
% % % % eventdata  structure with the following fields (see UIBUTTONGROUP)
% % % %	EventName: string 'SelectionChanged' (read only)
% % % %	OldValue: handle of the previously selected object or empty if none was selected
% % % %	NewValue: handle of the currently selected object
% % % % handles    structure with handles and user data (see GUIDATA)
% % % 
% % % switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
% % %     case 'radiobutton2a'
% % %         side_option = 'single';
% % %     case 'radiobutton2b'
% % %         side_option = 'double';
% % % end
% % % handles.metricdata.side_option = side_option;
% % % guidata(hObject,handles);


% % % % % --- Executes during object creation, after setting all properties.
% % % % function uipanel3_smooth_option_CreateFcn(hObject, eventdata, handles)
% % % % % hObject    handle to uipanel3_smooth_option (see GCBO)
% % % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % % handles    empty - handles not created until after all CreateFcns called
% % % % 
% % % % smooth_option = 'none';
% % % % handles.metricdata.smooth_option = smooth_option;
% % % % guidata(hObject,handles);
% % % % 
% % % % 
% % % % 
% % % % % --- Executes when selected object is changed in uipanel3_smooth_option.
% % % % function uipanel3_smooth_option_SelectionChangeFcn(hObject, eventdata, handles)
% % % % % hObject    handle to the selected object in uipanel3_smooth_option 
% % % % % eventdata  structure with the following fields (see UIBUTTONGROUP)
% % % % %	EventName: string 'SelectionChanged' (read only)
% % % % %	OldValue: handle of the previously selected object or empty if none was selected
% % % % %	NewValue: handle of the currently selected object
% % % % % handles    structure with handles and user data (see GUIDATA)
% % % % 
% % % % switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
% % % %     case 'radiobutton3a'
% % % %         smooth_option = 'none';
% % % %     case 'radiobutton3b'
% % % %         smooth_option = 'regular';
% % % %     case 'radiobutton3c'
% % % %         smooth_option = 'konnoOhmachi';
% % % % end
% % % % handles.metricdata.smooth_option = smooth_option;
% % % % guidata(hObject,handles);


% --- Executes on button press in pushbutton3_plot_all.
function pushbutton3_plot_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_plot_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.metricdata.step4_complete == 0
    msgbox('You haven''t selected any motions yet.','Warning');
else

    motion = handles.metricdata.motion;
    nr_motion = handles.metricdata.nr_motion;
    
    x_axis_scale = handles.metricdata.x_axis_scale;
    damping_ratio = handles.metricdata.damping_ratio;
    T_max = handles.metricdata.max_period;
    nr_pts = handles.metricdata.array_density;
    accel_unit = handles.metricdata.accel_unit;
    plot_option = 0;
    
    for i = 1 : 1 : nr_motion
        current_motion_filename = handles.metricdata.motion_file_name{i};
        current_motion = motion{i};
        
        [Tn,SA,PSA,SV,PSV,SD,fn] = responseSpectraExact(current_motion,...
                                  damping_ratio,T_max,plot_option,nr_pts);
        plotResponseSpectra(Tn,SA,PSA,SV,PSV,SD,current_motion,...
                            x_axis_scale,accel_unit,current_motion_filename);
        
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
    
    x_axis_scale = handles.metricdata.x_axis_scale;
    damping_ratio = handles.metricdata.damping_ratio;
    T_max = handles.metricdata.max_period;
    nr_pts = handles.metricdata.array_density;
    accel_unit = handles.metricdata.accel_unit;
    plot_option = 0;
    
    for i = 1 : 1 : nr_motion
        current_motion_filename = motion_filenames{i};
        current_motion = importdata(fullfile(motion_dir,current_motion_filename));
        
        [Tn,SA,PSA,SV,PSV,SD,fn] = responseSpectraExact(current_motion,...
                                  damping_ratio,T_max,plot_option,nr_pts);
        plotResponseSpectra(Tn,SA,PSA,SV,PSV,SD,current_motion,...
                            x_axis_scale,accel_unit,current_motion_filename);
    end
end


% --- Executes on button press in pushbutton5_calculate_all_and_save.
function pushbutton5_calculate_all_and_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_calculate_all_and_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.metricdata.step4_complete == 0
    msgbox('You haven''t selected any motions yet.','Warning');
else

    motion = handles.metricdata.motion;
    nr_motion = handles.metricdata.nr_motion;
    
    x_axis_scale = handles.metricdata.x_axis_scale;
    damping_ratio = handles.metricdata.damping_ratio;
    T_max = handles.metricdata.max_period;
    nr_pts = handles.metricdata.array_density;
    accel_unit = handles.metricdata.accel_unit;
    plot_option = 0;
    
    motion_dir = handles.metricdata.motion_dir_name;

    for i = 1 : 1 : nr_motion
        current_motion_filename = handles.metricdata.motion_file_name{i};
        current_motion = motion{i};
        
        [Tn,SA,PSA,SV,PSV,SD,fn] = responseSpectraExact(current_motion,...
                                  damping_ratio,T_max,plot_option,nr_pts);
        plotResponseSpectra(Tn,SA,PSA,SV,PSV,SD,current_motion,...
                            x_axis_scale,accel_unit,current_motion_filename);

        [str_tmp,motion_fname_without_ext,ext] = fileparts(current_motion_filename);
        RS_filename = sprintf('%s_response_spectra%s',motion_fname_without_ext,ext);
        dlmwrite(fullfile(motion_dir,RS_filename),[Tn,SA,PSA,SV,PSV,SD],'delimiter','\t','precision',6);
        
    end
    openFolder(motion_dir);
end


% --- Executes on button press in pushbutton6_calculate_selected_and_save.
function pushbutton6_calculate_selected_and_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6_calculate_selected_and_save (see GCBO)
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

    x_axis_scale = handles.metricdata.x_axis_scale;
    damping_ratio = handles.metricdata.damping_ratio;
    T_max = handles.metricdata.max_period;
    nr_pts = handles.metricdata.array_density;
    accel_unit = handles.metricdata.accel_unit;
    plot_option = 0;
    
    motion_dir = handles.metricdata.motion_dir_name;
    
    for i = 1 : 1 : nr_motion
        current_motion_filename = motion_filenames{i};
        current_motion = importdata(fullfile(motion_dir,current_motion_filename));
        
        [Tn,SA,PSA,SV,PSV,SD,fn] = responseSpectraExact(current_motion,...
                                  damping_ratio,T_max,plot_option,nr_pts);
        plotResponseSpectra(Tn,SA,PSA,SV,PSV,SD,current_motion,...
                            x_axis_scale,accel_unit,current_motion_filename);

        [str_tmp,motion_fname_without_ext,ext] = fileparts(current_motion_filename);
        RS_filename = sprintf('%s_response_spectra%s',motion_fname_without_ext,ext);
        dlmwrite(fullfile(motion_dir,RS_filename),[Tn,SA,PSA,SV,PSV,SD],'delimiter','\t','precision',6);
    end
    openFolder(motion_dir);
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;
