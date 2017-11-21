function varargout = SeismoSoil_Tools_PEER_To_2COL(varargin)
% SEISMOSOIL_TOOLS_PEER_TO_2COL M-file for SeismoSoil_Tools_PEER_To_2COL.fig
%      SEISMOSOIL_TOOLS_PEER_TO_2COL, by itself, creates a new SEISMOSOIL_TOOLS_PEER_TO_2COL or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_TOOLS_PEER_TO_2COL returns the handle to a new SEISMOSOIL_TOOLS_PEER_TO_2COL or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_TOOLS_PEER_TO_2COL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_TOOLS_PEER_TO_2COL.M with the given input arguments.
%
%      SEISMOSOIL_TOOLS_PEER_TO_2COL('Property','Value',...) creates a new SEISMOSOIL_TOOLS_PEER_TO_2COL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Tools_PEER_To_2COL_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Tools_PEER_To_2COL_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Tools_PEER_To_2COL

% Last Modified by GUIDE v2.5 21-May-2015 15:16:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Tools_PEER_To_2COL_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Tools_PEER_To_2COL_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Tools_PEER_To_2COL is made visible.
function SeismoSoil_Tools_PEER_To_2COL_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Tools_PEER_To_2COL (see VARARGIN)

% Choose default command line output for SeismoSoil_Tools_PEER_To_2COL
handles.output = hObject;

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Tools_PEER_To_2COL wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Tools_PEER_To_2COL_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


% --- Executes during object creation, after setting all properties.
function pushbutton2_select_motions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2_select_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

try
    handles.metricdata.uigetfile_start_dir;
catch
    handles.metricdata.uigetfile_start_dir = pwd;
end

handles.metricdata.step4_complete = 0;
guidata(hObject,handles);



% --- Executes on button press in pushbutton2_select_motions.
function pushbutton2_select_motions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_select_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;

global start_dir0;

filter_spec = {'*.AT2','NGA Acceleration Format (*.AT2)';'*.*','All Files (*.*)'};
dlg_title = 'Select input file(s)...';
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
    motion{i} = importdata(fullfile(motion_dir_name,motion_file_name{i}),...
                           ' ',4);
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


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


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


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

% --- Executes during object creation, after setting all properties.
function uipanel1_unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.factor_to_SI = 1;
guidata(hObject,handles);





% --- Executes during object creation, after setting all properties.
function uipanel7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.factor_from_SI = 1;
guidata(hObject,handles);




% --- Executes when selected object is changed in uipanel7.
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel7 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton5a'
        factor_from_SI = 1;
    case 'radiobutton5b'
        factor_from_SI = 100;
    case 'radiobutton5c'
        factor_from_SI = 1/9.81;
end
handles.metricdata.factor_from_SI = factor_from_SI;
guidata(hObject,handles);




% --- Executes on button press in pushbutton3_convert_all.
function pushbutton3_convert_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_convert_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;

if handles.metricdata.step4_complete == 0
    msgbox('You haven''t selected any files yet.','Warning');
else
   
    motion = handles.metricdata.motion;
    nr_motion = handles.metricdata.nr_motion;
    
    motion_dir_name = handles.metricdata.motion_dir_name;
    for i = 1 : 1 : nr_motion
        current_motion_filename = handles.metricdata.motion_file_name{i};
        [~,fname,ext] = fileparts(current_motion_filename);
        fprintf('%s\n',current_motion_filename);
        
        switch handles.metricdata.factor_from_SI % Get Tag of selected object.
            case 1
                str1 = 'SI';
            case 100
                str1 = 'gal';
            case 1/9.81
                str1 = 'g';
        end
        new_fname = sprintf('%s_(unit=%s).txt',fname,str1);
        
        current_motion_struct = motion{i};
        factor_from_SI = handles.metricdata.factor_from_SI;
        factor_to_SI = 9.81;
        
        accel_matrix = current_motion_struct.data;
        
        if any(strcmp('colheaders',fields(current_motion_struct)))
            dt = str2double(current_motion_struct.colheaders{4});
        else  % if current_motion_struct does not have a field named 'colheaders'
            dt_info = current_motion_struct.textdata{end};
            dt_str = dt_info(18:26);
            dt = str2double(dt_str);
        end
        
        accel_matrix_tr = transpose(accel_matrix);
        npts = numel(accel_matrix);
        
        time = (dt : dt : dt*npts)';
        accel = zeros(npts,1);
        
        for j = 1 : 1 : npts
            accel(j) = accel_matrix_tr(j);
        end
        
        accel = accel * factor_to_SI * factor_from_SI;
        
        dlmwrite(fullfile(motion_dir_name,new_fname),[time,accel],'delimiter','\t','precision',6);
        
    end
    
    choice = questdlg('Finished. Open containing folder?', ...
	'Finished', ...
	'Yes','No','No');
    switch choice
        case 'Yes'
            dir_absolute = cd(cd(motion_dir_name));
            command_text = sprintf('explorer.exe %s',dir_absolute);
            system(command_text);
        case 'No'
            % do nothing
    end
    
end

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

% --- Executes on button press in pushbutton4_convert_selected.
function pushbutton4_convert_selected_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_convert_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;

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

    motion_dir_name = handles.metricdata.motion_dir_name;
    for i = 1 : 1 : nr_motion
        current_motion_filename = motion_filenames{i};
        [~,fname,ext] = fileparts(current_motion_filename);
        
        switch handles.metricdata.factor_from_SI % Get Tag of selected object.
            case 1
                str1 = 'SI';
            case 100
                str1 = 'gal';
            case 1/9.81
                str1 = 'g';
        end
        new_fname = sprintf('%s_in_%s%s',fname,str1,ext);
        
        current_motion = importdata(fullfile(motion_dir_name,current_motion_filename));
        factor_to_SI = handles.metricdata.factor_to_SI;
        factor_from_SI = handles.metricdata.factor_from_SI;
        
        time = current_motion(:,1);
        accel = current_motion(:,2);
        
        accel = accel * factor_to_SI * factor_from_SI;
        
        dlmwrite(fullfile(motion_dir_name,new_fname),[time,accel],'delimiter','\t','precision',6);
        
    end
    
    choice = questdlg('Finished. Open containing folder?', ...
        'Finished', ...
        'Yes','No','No');
    switch choice
        case 'Yes'
            dir_absolute = cd(cd(motion_dir_name));
            command_text = sprintf('explorer.exe %s',dir_absolute);
            system(command_text);
        case 'No'
            % do nothing
    end
end


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


% --- Executes on button press in pushbutton12_close_all.
function pushbutton12_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;



% --- Executes on button press in pushbutton1_return_to_tools.
function pushbutton1_return_to_tools_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_return_to_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;
close SeismoSoil_Tools_PEER_To_2COL;
SeismoSoil_Tools;


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.metricdata.dt_entered = 0;
guidata(hObject,handles);
