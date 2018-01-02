function varargout = SeismoSoil_Tools_Linear_Ampl_Factor(varargin)
% SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR M-file for SeismoSoil_Tools_Linear_Ampl_Factor.fig
%      SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR, by itself, creates a new SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR returns the handle to a new SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR.M with the given input arguments.
%
%      SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR('Property','Value',...) creates a new SEISMOSOIL_TOOLS_LINEAR_AMPL_FACTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Tools_Linear_Ampl_Factor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Tools_Linear_Ampl_Factor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Tools_Linear_Ampl_Factor

% Last Modified by GUIDE v2.5 01-Jan-2018 18:09:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Tools_Linear_Ampl_Factor_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Tools_Linear_Ampl_Factor_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Tools_Linear_Ampl_Factor is made visible.
function SeismoSoil_Tools_Linear_Ampl_Factor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Tools_Linear_Ampl_Factor (see VARARGIN)

% Choose default command line output for SeismoSoil_Tools_Linear_Ampl_Factor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);

% UIWAIT makes SeismoSoil_Tools_Linear_Ampl_Factor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Tools_Linear_Ampl_Factor_OutputFcn(hObject, eventdata, handles) 
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

close SeismoSoil_Tools_Linear_Ampl_Factor;
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
    checkInputs(motion{i},'vs_profile',sprintf('Vs profile #%d',i));
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
        [freq_array,AF_ro,TF_ro,f0_ro,AF_in,TF_in,AF_bh,TF_bh,f0_bh]...
            = tfLIN(current_profile,'off');
        plotProfileAndLAF(current_profile,freq_array,AF_ro,f0_ro,AF_bh,f0_bh,...
            current_profile_filename);
    end
    
end


function plotProfileAndLAF(vs_profile,freq,AF_ro,f0_ro,AF_bh,f0_bh,name)

hfig = figure;
set(hfig,'unit','inches','paperposition',[3,3,7.5,5.5]);

subplot(2,3,[1,4]);
[x,y] = genProfilePlotArrayFromMatrix(vs_profile);
plot(x,y,'k','linewidth',1.5);
set(gca,'Ydir','reverse','fontsize',10);
%ylim([0, sum(vs_profile(:,1))]);
ylim([0, inf]);
xlim([0, max(vs_profile(:,2))*1.2]);
xlabel('Vs [m/s]','fontsize',10);
ylabel('Depth [m]','fontsize',10);
grid on;
title(name,'interpreter','none');

subplot(2,3,2:3);
semilogx(freq,AF_ro,'k','linewidth',1.5); hold on;
plot(f0_ro,AF_ro(find(freq==f0_ro)),'ro','linewidth',1.75);
annotation('textbox',[.42 .8 .1 .1],'String',sprintf('f_0 = %.2f Hz',f0_ro),'fontsize',9,'backgroundcolor','w');
xlim([min(freq) max(freq)]);
ylim([0 max(AF_ro)*1.15]);
ylabel('Amplification Factor','fontsize',10);
grid on;
title('Rock outcrop');

subplot(2,3,5:6);
semilogx(freq,AF_bh,'k','linewidth',1.5); hold on;
plot(f0_bh,AF_bh(find(freq==f0_bh)),'ro','linewidth',1.75);
annotation('textbox',[.42 .32 .1 .1],'String',sprintf('f_0 = %.2f Hz',f0_bh),'fontsize',9,'backgroundcolor','w');
xlabel('Frequency [Hz]','fontsize',10);
xlim([min(freq) max(freq)]);
ylim([0 max(AF_bh)*1.15]);
ylabel('Amplification Factor','fontsize',10);
grid on;
title('Borehole');


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
        current_profile_filename = motion_filenames{i};
        current_profile = importdata(fullfile(motion_dir,current_profile_filename));
        [freq_array,AF_ro,TF_ro,f0_ro,AF_in,TF_in,AF_bh,TF_bh,f0_bh]...
            = tfLIN(current_profile,'off');
        plotProfileAndLAF(current_profile,freq_array,AF_ro,f0_ro,AF_bh,f0_bh,...
            current_profile_filename);
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
