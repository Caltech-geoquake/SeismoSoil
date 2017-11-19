function varargout = SeismoSoil_Tools_Motion_Plotter(varargin)
% SEISMOSOIL_TOOLS_MOTION_PLOTTER M-file for SeismoSoil_Tools_Motion_Plotter.fig
%      SEISMOSOIL_TOOLS_MOTION_PLOTTER, by itself, creates a new SEISMOSOIL_TOOLS_MOTION_PLOTTER or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_TOOLS_MOTION_PLOTTER returns the handle to a new SEISMOSOIL_TOOLS_MOTION_PLOTTER or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_TOOLS_MOTION_PLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_TOOLS_MOTION_PLOTTER.M with the given input arguments.
%
%      SEISMOSOIL_TOOLS_MOTION_PLOTTER('Property','Value',...) creates a new SEISMOSOIL_TOOLS_MOTION_PLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Tools_Motion_Plotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Tools_Motion_Plotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Tools_Motion_Plotter

% Last Modified by GUIDE v2.5 25-Jul-2017 23:19:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Tools_Motion_Plotter_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Tools_Motion_Plotter_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Tools_Motion_Plotter is made visible.
function SeismoSoil_Tools_Motion_Plotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Tools_Motion_Plotter (see VARARGIN)

% Choose default command line output for SeismoSoil_Tools_Motion_Plotter
handles.output = hObject;

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Tools_Motion_Plotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Tools_Motion_Plotter_OutputFcn(hObject, eventdata, handles) 
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
    
    [check_flag,err_msg] = checkInputs(motion{i},'motion');
    if check_flag == -1
        if nr_motion > 1  % if user loads more than one ground motions
            err_msg = sprintf('Motion #%d %s',i,err_msg(6:end));
        end
        fprintf('***** %s *****\n',err_msg);
        msgbox(err_msg, 'Warning');
    end
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

accel_unit_division_factor = 1;  % default value of unit selection
input_accel_unit_ = 'm/s^2'; % default value of unit selection
handles.metricdata.accel_unit_division_factor = accel_unit_division_factor;
handles.metricdata.input_accel_unit = input_accel_unit_;
guidata(hObject,handles);



% --- Executes when selected object is changed in uipanel1_unit.
function uipanel1_unit_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1_unit 
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

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


% --- Executes during object creation, after setting all properties.
function checkbox1_arias_intensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox1_arias_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

arias_intensity_option = 0;
handles.metricdata.arias_intensity_option = arias_intensity_option;
guidata(hObject,handles);




% --- Executes on button press in checkbox1_arias_intensity.
function checkbox1_arias_intensity_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1_arias_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1_arias_intensity

arias_intensity_option = get(hObject,'value');
handles.metricdata.arias_intensity_option = arias_intensity_option;
guidata(hObject,handles);



% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


% --- Executes during object creation, after setting all properties.
function checkbox2_show_RMS_accel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox2_show_RMS_accel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

show_RMS_accel = 0;
handles.metricdata.show_RMS_accel = show_RMS_accel;
guidata(hObject,handles);



% --- Executes on button press in checkbox2_show_RMS_accel.
function checkbox2_show_RMS_accel_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2_show_RMS_accel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2_show_RMS_accel

show_RMS_accel = get(hObject,'value');
handles.metricdata.show_RMS_accel = show_RMS_accel;
guidata(hObject,handles);




% --- Executes on button press in pushbutton13_plot_PGA_histogram.
function pushbutton13_plot_PGA_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13_plot_PGA_histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.metricdata.step4_complete == 0
    msgbox('You haven''t selected any motions yet.','Warning');
else

    motion = handles.metricdata.motion;
    nr_motion = handles.metricdata.nr_motion;
    unit = handles.metricdata.input_accel_unit;
    
    pga = zeros(nr_motion,1);
    for i = 1 : 1 : nr_motion
        current_motion = motion{i};    
        pga(i) = max(abs(current_motion(:,2)));
    end
    
    if nr_motion <= 20
        nr_bins = 10;
    elseif nr_motion <= 50
        nr_bins = 20;
    elseif nr_motion <= 100
        nr_bins = 30;
    else
        nr_bins = 40;
    end
    
    figure;
    histogram(pga,nr_bins);
    xlabel(sprintf('PGA [%s]',unit));
    ylabel('Number of ground motions');
    grid on;
    title(sprintf('Total: %d motions',nr_motion));
    
end


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

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

    for i = 1 : 1 : nr_motion
        current_motion_filename = handles.metricdata.motion_file_name{i};
        current_motion = motion{i};
        unit = handles.metricdata.input_accel_unit;
        factor = handles.metricdata.accel_unit_division_factor;
        show_RMS_accel = handles.metricdata.show_RMS_accel;
        plotMotion(current_motion,unit,factor,current_motion_filename,show_RMS_accel);
        
        if handles.metricdata.arias_intensity_option == 1
            motion_in_SI = current_motion;
            motion_in_SI(:,2) = motion_in_SI(:,2)/factor;
            Ia = calcAriasIntensity(motion_in_SI);
            fh = figure;
            width = 5; height = 3; 
            set(fh,'units','inches','position',[5.1+(10-width)/2,4+(6-height)/2,width,height]);
            plot(Ia(:,1),Ia(:,2),'linewidth',1.5);
            xlabel('Time (s)','fontsize',12);
            ylabel('Arias Intensity (m/s)','fontsize',12);
            xlim([0 max(Ia(:,1))]);
            grid on;
        end
    end
    
end

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

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
        unit = handles.metricdata.input_accel_unit;
        factor = handles.metricdata.accel_unit_division_factor;
        show_RMS_accel = handles.metricdata.show_RMS_accel;
        plotMotion(current_motion,unit,factor,current_motion_filename,show_RMS_accel);
        
        if handles.metricdata.arias_intensity_option == 1
            motion_in_SI = current_motion;
            motion_in_SI(:,2) = motion_in_SI(:,2)/factor;
            Ia = calcAriasIntensity(motion_in_SI);
            fh = figure;
            width = 5; height = 3; 
            set(fh,'units','inches','position',[5.1+(10-width)/2,4+(6-height)/2,width,height]);
            plot(Ia(:,1),Ia(:,2),'linewidth',1.5);
            xlabel('Time (s)','fontsize',12);
            ylabel('Arias Intensity (m/s)','fontsize',12);
            xlim([0 max(Ia(:,1))]);
            grid on;
        end
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

close SeismoSoil_Tools_Motion_Plotter;
SeismoSoil_Tools;

