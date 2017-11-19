function varargout = SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis(varargin)
% SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS MATLAB code for SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis.fig
%      SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS, by itself, creates a new SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS returns the handle to a new SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS.M with the given input arguments.
%
%      SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS('Property','Value',...) creates a new SEISMOSOIL_EQUIVALENT_LINEAR_FREQUENCY_DEPENDENT_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis

% Last Modified by GUIDE v2.5 24-Aug-2014 23:38:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis is made visible.
function SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis (see VARARGIN)

clc;

% Choose default command line output for SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1)

% UIWAIT makes SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function pushbutton1a_select_profile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1a_select_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.step1_complete = 0;
guidata(hObject,handles);


% --- Executes on button press in pushbutton1a_select_profile.
function pushbutton1a_select_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1a_select_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select Vs profile data file...';
[profile_file_name,profile_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');

if isequal(profile_file_name,0)
    disp('User did not select any files.');
else
    disp(['User selected ',fullfile(profile_dir_name,profile_file_name)]);
    start_dir0 = profile_dir_name;

    handles.metricdata.profile_file_name = profile_file_name;
    handles.metricdata.profile_dir_name = profile_dir_name;
    
    % rename profile
    if findstr('profile_',profile_file_name) % if the profile filename contains "profile_"
        sitecode_ = profile_file_name(9:end-4);
    else
        sitecode_ = profile_file_name;
    end
    handles.metricdata.sitecode = sitecode_;
    set(handles.edit1a_profile_name,'string',sitecode_);
    
    profile_ = importdata(fullfile(profile_dir_name,profile_file_name));
    handles.metricdata.profile = profile_;
    
    handles.metricdata.step1_complete = 1;
    
    folder_name = sprintf('%s_Freq_Dependent_Equiv_Linear_Results',sitecode_);
    output_dir_ = fullfile('./',folder_name);
    handles.metricdata.output_dir = output_dir_;
end

guidata(hObject,handles);



% --- Executes on selection change in popupmenu1a_bedrock_type.
function popupmenu1a_bedrock_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1a_bedrock_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1a_bedrock_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1a_bedrock_type
contents = cellstr(get(hObject,'String'));
bedrock_type_ = contents{get(hObject,'Value')};
handles.metricdata.bedrock_type = bedrock_type_;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1a_bedrock_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1a_bedrock_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

contents = cellstr(get(hObject,'String'));
bedrock_type_ = contents{get(hObject,'Value')};
handles.metricdata.bedrock_type = bedrock_type_;
guidata(hObject,handles);


% --- Executes on button press in pushbutton1b_plot_profile.
function pushbutton1b_plot_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1b_plot_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

profile_file_name_ = handles.metricdata.profile_file_name;
profile_dir_name_ = handles.metricdata.profile_dir_name;

full_filename = fullfile(profile_dir_name_,profile_file_name_);
plotVsProfileFromFilename(full_filename);
title(profile_file_name_,'interpreter','none');


% --- Executes during object creation, after setting all properties.
function edit1a_profile_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1a_profile_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit1a_profile_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit1a_profile_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1a_profile_name as text
%        str2double(get(hObject,'String')) returns contents of edit1a_profile_name as a double


sitecode_ = get(hObject,'String');
handles.metricdata.sitecode = sitecode_;

folder_name = sprintf('%s_Freq_Dependent_Equiv_Linear_Results',sitecode_);
output_dir_ = fullfile('./',folder_name);
handles.metricdata.output_dir = output_dir_;

guidata(hObject,handles);


% --- Executes on button press in pushbutton1c_help_dlgbox.
function pushbutton1c_help_dlgbox_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1c_help_dlgbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox({'If you don''t specify a Vs profile name, a site code modified',...
        'from the file name is used. For example, the file name',...
        'is usually "profile_CE.11023.dat", then the Vs profile name',...
        'will be "CE.11023".'},...
        'Help');


% --- Executes during object creation, after setting all properties.
function pushbutton2a_select_curve_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2a_select_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.metricdata.step2_complete = 0;
guidata(hObject,handles);
    
    
% --- Executes on button press in pushbutton2a_select_curve.
function pushbutton2a_select_curve_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2a_select_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*', 'All Files (*.*)'};
dlg_title = 'Select "curve" data file...';
[curve_file_name,curve_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');

if isequal(curve_file_name,0)
    disp('User did not select any files.');
else
    disp(['User selected ',fullfile(curve_dir_name,curve_file_name)]);
    start_dir0 = curve_dir_name; % update "start dir"

    handles.metricdata.curve_file_name = curve_file_name;
    handles.metricdata.curve_dir_name = curve_dir_name;
    
    curve = importdata(fullfile(curve_dir_name,curve_file_name));
    handles.metricdata.curve = curve;
    
    handles.metricdata.step2_complete = 1;
end

guidata(hObject,handles);


% --- Executes on button press in pushbutton2b_plot_curves.
function pushbutton2b_plot_curves_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2b_plot_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curve_matrix = handles.metricdata.curve;
plotCurves(curve_matrix);
title(handles.metricdata.curve_file_name,'interpreter','none');



% --- Executes during object creation, after setting all properties.
function uipanel13_rho_unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel13_rho_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

rho_unit_division_factor = 1;  % default value of rho unit factor
handles.metricdata.rho_unit_division_factor = rho_unit_division_factor;
guidata(hObject,handles);



% --- Executes when selected object is changed in uipanel13_rho_unit.
function uipanel13_rho_unit_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel13_rho_unit 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1a'
        rho_unit_division_factor = 1;
    case 'radiobutton1b'
        rho_unit_division_factor = 1/1000;
end
handles.metricdata.rho_unit_division_factor = rho_unit_division_factor;
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function uipanel1b_xi_unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1b_xi_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

xi_unit_division_factor = 1;  % default value
handles.metricdata.xi_unit_division_factor = xi_unit_division_factor;
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel1b_xi_unit.
function uipanel1b_xi_unit_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1b_xi_unit 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1c'
        xi_unit_division_factor = 1;
    case 'radiobutton1d'
        xi_unit_division_factor = 100;
end
handles.metricdata.xi_unit_division_factor = xi_unit_division_factor;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function uipanel_unit_accel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_unit_accel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

accel_unit_division_factor = 1;  % default value of unit selection
input_accel_unit_ = 'm/s^2'; % default value of unit selection
handles.metricdata.accel_unit_division_factor = accel_unit_division_factor;
handles.metricdata.input_accel_unit = input_accel_unit_;
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel_unit_accel.
function uipanel_unit_accel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_unit_accel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton4a'
        accel_unit_division_factor = 1;
        input_accel_unit_ = 'm/s^2';
    case 'radiobutton4b'
        accel_unit_division_factor = 100;
        input_accel_unit_ = 'gal';
    case 'radiobutton4c'
        accel_unit_division_factor = 1/9.81;
        input_accel_unit_ = 'g';
end
handles.metricdata.accel_unit_division_factor = accel_unit_division_factor;
handles.metricdata.input_accel_unit = input_accel_unit_;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function pushbutton4a_select_motion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton4a_select_motion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.step3_complete = 0;
guidata(hObject,handles);


% --- Executes on button press in pushbutton4a_select_motion.
function pushbutton4a_select_motion_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4a_select_motion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*', 'All Files (*.*)'};
dlg_title = 'Select input motion(s)...';
[motion_file_name,motion_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','on');

if isequal(motion_file_name,0)
    disp('User did not select any files.');
else
    if ischar(motion_file_name) % if motion_file_name is a string
        disp('User selected 1 motion.');
    else
        disp(['User selected ',mat2str(length(motion_file_name)),' motions.']);
    end
    start_dir0 = motion_dir_name; % update "start dir"
    
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
    handles.metricdata.step3_complete = 1;
    
    % Initialize baseline correction log
    baseline_correction_log = zeros(nr_motion,1);
    handles.metricdata.baseline_correction_log = baseline_correction_log;
    
    % Initialize listbox_motions and set selected_motion_indices to 1
    handles.metricdata.selected_motion_indices = 1;
    temp = handles.metricdata.motion_file_name;
    set(handles.listbox_motions,'String',temp,'Value',1);
    handles.metricdata.motion_listbox_contents = temp;
    
end

guidata(hObject,handles);


% --- Executes on button press in pushbutton4b_plot_all_motions.
function pushbutton4b_plot_all_motions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4b_plot_all_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.metricdata.step3_complete == 0
    msgbox('You haven''t selected any motions yet.','Warning');
else

    motion = handles.metricdata.motion;
    nr_motion = handles.metricdata.nr_motion;

    for i = 1 : 1 : nr_motion
        current_motion_filename = handles.metricdata.motion_file_name{i};
        current_motion = motion{i};
        unit = handles.metricdata.input_accel_unit;
        factor = handles.metricdata.accel_unit_division_factor;
        plotMotion(current_motion,unit,factor,current_motion_filename);
    end
    
end


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


% --- Executes on selection change in listbox_motions.
function listbox_motions_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_motions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_motions

motion_listbox_contents = get(hObject,'String');  
handles.metricdata.motion_listbox_contents = motion_listbox_contents;
% motion_selected = motion_listbox_contents{get(hObject,'Value')};
% % selected_motion_index = get(hObject,'value');
% % if isfield(handles.metricdata,'selected_motion_index_array')
% %     selected_motion_index_array = handles.metricdata.selected_motion_index_array;
% %     selected_motion_index_array = [selected_motion_index_array;{selected_motion_index}];
% %     % update 'selected_motion_index_array'
% %     handles.metricdata.selected_motion_index_array = selected_motion_index_array;
% % else
% %     handles.metricdata.selected_motion_index_array = selected_motion_index;
% % end
selected_motion_indices = get(hObject,'value');
handles.metricdata.selected_motion_indices = selected_motion_indices;
guidata(hObject,handles);


% --- Executes on button press in pushbutton4c_plot_selected_motions.
function pushbutton4c_plot_selected_motions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4c_plot_selected_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.metricdata.step3_complete == 0
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
        unit = handles.metricdata.input_accel_unit;
        factor = handles.metricdata.accel_unit_division_factor;
        plotMotion(current_motion,unit,factor,current_motion_filename);
    end
end


% --- Executes on button press in pushbutton4e_baseline_correction_selected.
function pushbutton4e_baseline_correction_selected_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4e_baseline_correction_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cutoff_freq = 0.2;
baseline_correction_log = handles.metricdata.baseline_correction_log;

if handles.metricdata.nr_motion == 1
    selected_motion_indices = 1;
else
    selected_motion_indices = handles.metricdata.selected_motion_indices;
end

flag = 0;
if prod(baseline_correction_log(selected_motion_indices)) == 1
    msgbox({'All of the selected motions have been corrected.',...
            'They will not be processed again.'},'Information');
    flag = -1;
else if sum(baseline_correction_log(selected_motion_indices)) > 0 % there are 1's in the array
    msgbox({'Some of the selected motions have been corrected.',...
            'They will not be processed again.'},'Information');
    end
end

for i = selected_motion_indices
    if baseline_correction_log(i) == 0
        accel = handles.metricdata.motion{i};
        baseline_result = baselineCorrHighPass(accel,cutoff_freq,1);
        handles.metricdata.motion{i} = baseline_result;
        baseline_correction_log(i) = 1; % marks that the i-th motion has been corrected
    end
end
if flag == -1
else
    msgbox('All baseline correction finished.','Finished');
end
 
handles.metricdata.baseline_correction_log = baseline_correction_log;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton4d_baseline_correction_all_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton4d_baseline_correction_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%
% handles.metricdata.hit_baseline_button_counter = 0;
% % % This counter is to be used when this button is pushed. If it has been hit
% % % befored, then re-hitting this button takes no effect. This is to ensure
% % % that the user does not "over-correct" the motion.
% guidata(hObject,handles);

% --- Executes on button press in pushbutton4d_baseline_correction_all.
function pushbutton4d_baseline_correction_all_Callback(hObject,eventdata,handles)
% hObject    handle to pushbutton4d_baseline_correction_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flag = 0;
cutoff_freq = 0.2;
nr_motion = handles.metricdata.nr_motion;
baseline_correction_log = handles.metricdata.baseline_correction_log;
if sum(baseline_correction_log) == 0 % if all elements is zero
    for i = 1 : 1 : nr_motion
        accel = handles.metricdata.motion{i};
        baseline_result = baselineCorrHighPass(accel,cutoff_freq,1);
        handles.metricdata.motion{i} = baseline_result;
        baseline_correction_log(i) = 1; % marks that the i-th motion has been corrected
    end
    handles.metricdata.baseline_correction_log = baseline_correction_log;
    guidata(hObject,handles);
else if prod(baseline_correction_log) == 1 % they've all been processed
        msgbox({'All of the motions have already been corrected.',...
                'They will not be processed again.'},'Information');
        flag = -1;
    else % one or some of the motions has already been corrected
        msgbox({'Some of the motions you selected have already been corrected.',...
                'They will not be processed again.'},'Information');
        for i = 1 : 1 : nr_motion
            if baseline_correction_log(i) == 0
                accel = handles.metricdata.motion{i};
                baseline_result = baselineCorrHighPass(accel,cutoff_freq,1);
                handles.metricdata.motion{i} = baseline_result;
                baseline_correction_log(i) = 1; % marks that the i-th motion has been corrected
            end
        end
        handles.metricdata.baseline_correction_log = baseline_correction_log;
        guidata(hObject,handles);
    end
end
if flag == -1
else
    msgbox('All baseline correction finished.','Finished');
end

% --- Executes during object creation, after setting all properties.
function uipanel4b_motion_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel4b_motion_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

motion_type_ = 'incident'; % default value
handles.metricdata.motion_type = motion_type_;
guidata(hObject,handles);

% --- Executes when selected object is changed in uipanel4b_motion_type.
function uipanel4b_motion_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4b_motion_type 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton441'
        motion_type_ = 'incident';
    case 'radiobutton442'
        motion_type_ = 'borehole';
    case 'radiobutton443'
        motion_type_ = 'outcrop';
end
handles.metricdata.motion_type = motion_type_;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton5a_select_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton5a_select_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% % sitecode_ = handles.metricdata.sitecode;
% % folder_name = sprintf('%s_NL_result',sitecode_);
% % output_dir_ = fullfile('./',folder_name);
% % mkdir(output_dir);
% % handles.metricdata.output_dir = output_dir_;
% % guidata(hObject,handles);



% --- Executes on button press in pushbutton5a_select_folder.
function pushbutton5a_select_folder_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5a_select_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;
output_dir_ = uigetdir(start_dir0,'Select directory for exporting results...');
fprintf('Output directory: %s\n',output_dir_);
handles.metricdata.output_dir = output_dir_;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function checkbox1_open_result_folder_option_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox1_open_result_folder_option (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.open_result_folder_option = 1; % default: checked
guidata(hObject,handles);



% --- Executes on button press in checkbox1_open_result_folder_option.
function checkbox1_open_result_folder_option_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1_open_result_folder_option (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1_open_result_folder_option

open_result_folder_option = get(hObject,'Value');
handles.metricdata.open_result_folder_option = open_result_folder_option;
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function checkbox2_view_results_as_popout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox2_view_results_as_popout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.view_results_as_popup_option = 0; % default: unchecked
guidata(hObject,handles);




% --- Executes on button press in checkbox2_view_results_as_popout.
function checkbox2_view_results_as_popout_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2_view_results_as_popout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2_view_results_as_popout

view_results_as_popup_option = get(hObject,'Value');
handles.metricdata.view_results_as_popup_option = view_results_as_popup_option;
guidata(hObject,handles);



% --- Executes on button press in pushbutton_start_NL.
function pushbutton_start_NL_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start_NL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gcp;  % start parallel workers

warning('off','all');

step1 = handles.metricdata.step1_complete;
step2 = handles.metricdata.step2_complete;
step3 = handles.metricdata.step3_complete;
% step4 = handles.metricdata.step4_complete;
% step = [step1,step2,step3,step4];
step = [step1,step2,step3];

str1 = ' ';
if step1*step2*step3 == 0 % if a certain step is not completed
    for i = 1 : 1 : 3
        if step(i) == 0
            str1 = [str1,num2str(i),','];
        end
    end
    str1 = str1(1:end-1); % get rid of the last comma
    warning_text = sprintf('You haven''t finished step(s)%s.\nYou need to finish it/them before\nstarting the analysis.',str1);
    h_msgbox = msgbox(warning_text,'Warning');
else
    vs_profile = handles.metricdata.profile;
    curve = handles.metricdata.curve;
%     H2n = handles.metricdata.H2n;
    nr_motion = handles.metricdata.nr_motion;
    motion = handles.metricdata.motion;
    motion_name = handles.metricdata.motion_file_name;
    output_dir = handles.metricdata.output_dir;
    factor_rho = handles.metricdata.rho_unit_division_factor;
    factor_xi = handles.metricdata.xi_unit_division_factor;
    unit_factor_accel = handles.metricdata.accel_unit_division_factor;
    bedrock_type = handles.metricdata.bedrock_type;
    motion_type = handles.metricdata.motion_type;
    
%     mkdir(output_dir);
    
    if handles.metricdata.view_results_as_popup_option == 1
        fig_visible_option = 'on';
    else
        fig_visible_option = 'off';
    end
    
    tic;
    [ok_to_proceed,h_running] = runEquivLinearFreqDepFromGUI(vs_profile,...
        curve,nr_motion,motion,motion_name,output_dir,...
        factor_rho,factor_xi,unit_factor_accel,bedrock_type,...
        motion_type,fig_visible_option);
    elapsed_time = toc;
    fprintf('Simulation finished. Elapsed time: %.1f sec.\n',elapsed_time);
    
end

if ok_to_proceed == 1
    if ishandle(h_running) % if the users haven't close the figure
        close(h_running);
    end
    msgbox(sprintf('Analysis complete!\nElapsed time: %1f sec.',elapsed_time),'Finished');
    if handles.metricdata.open_result_folder_option == 1
        output_dir_absolute = cd(cd(output_dir));
        if ismac()
            output_dir_absolute = regexprep(output_dir_absolute,' ','\\ ');
        end
        
        command_Windows = sprintf('explorer.exe %s',output_dir_absolute);
        command_Mac = sprintf('open %s',output_dir_absolute);
        if ispc()
            system(command_Windows);
        elseif ismac()
            system(command_Mac);
        end
    end
end


% --- Executes on button press in pushbutton17_open_result_folder.
function pushbutton17_open_result_folder_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17_open_result_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

output_dir = handles.metricdata.output_dir;
output_dir_absolute = cd(cd(output_dir));
if ismac()
    output_dir_absolute = regexprep(output_dir_absolute,' ','\\ ');
end

command_Windows = sprintf('explorer.exe %s',output_dir_absolute);
command_Mac = sprintf('open %s',output_dir_absolute);
if ispc()
    system(command_Windows);
elseif ismac()
    system(command_Mac);
end


% --- Executes on button press in pushbutton0_return_to_main.
function pushbutton0_return_to_main_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton0_return_to_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Equivalent_Linear_Frequency_Dependent_Analysis;
SeismoSoil;


% --- Executes on button press in pushbutton19_close_all.
function pushbutton19_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;
