function varargout = SeismoSoil_Tools_Deconvolution(varargin)
% SEISMOSOIL_TOOLS_DECONVOLUTION MATLAB code for SeismoSoil_Tools_Deconvolution.fig
%      SEISMOSOIL_TOOLS_DECONVOLUTION, by itself, creates a new SEISMOSOIL_TOOLS_DECONVOLUTION or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_TOOLS_DECONVOLUTION returns the handle to a new SEISMOSOIL_TOOLS_DECONVOLUTION or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_TOOLS_DECONVOLUTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_TOOLS_DECONVOLUTION.M with the given input arguments.
%
%      SEISMOSOIL_TOOLS_DECONVOLUTION('Property','Value',...) creates a new SEISMOSOIL_TOOLS_DECONVOLUTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Tools_Deconvolution_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Tools_Deconvolution_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Tools_Deconvolution

% Last Modified by GUIDE v2.5 07-Jun-2016 17:23:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Tools_Deconvolution_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Tools_Deconvolution_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Tools_Deconvolution is made visible.
function SeismoSoil_Tools_Deconvolution_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Tools_Deconvolution (see VARARGIN)

% Choose default command line output for SeismoSoil_Tools_Deconvolution
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);

% UIWAIT makes SeismoSoil_Tools_Deconvolution wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Tools_Deconvolution_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_return.
function pushbutton1_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Tools_Deconvolution;
SeismoSoil_Tools;



function edit1_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_thickness as text
%        str2double(get(hObject,'String')) returns contents of edit1_thickness as a double
handles.metricdata.thickness = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit1_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.metricdata.thickness = 0;
handles.metricdata.thickness_handle = hObject;

guidata(hObject,handles);


function edit2_Vs_Callback(hObject, eventdata, handles)
% hObject    handle to edit2_Vs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2_Vs as text
%        str2double(get(hObject,'String')) returns contents of edit2_Vs as a double
handles.metricdata.Vs = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit2_Vs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2_Vs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.metricdata.Vs = 0;
handles.metricdata.Vs_handle = hObject;

guidata(hObject,handles);


function edit3_density_Callback(hObject, eventdata, handles)
% hObject    handle to edit3_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3_density as text
%        str2double(get(hObject,'String')) returns contents of edit3_density as a double
handles.metricdata.density = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit3_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.metricdata.density = 0;
handles.metricdata.density_handle = hObject;

guidata(hObject,handles);


function edit4_damping_Callback(hObject, eventdata, handles)
% hObject    handle to edit4_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4_damping as text
%        str2double(get(hObject,'String')) returns contents of edit4_damping as a double
handles.metricdata.damping = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit4_damping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.metricdata.damping = 0;
handles.metricdata.damping_handle = hObject;

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton2_select_surface_motion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2_select_surface_motion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.select_motion_complete = 0;

guidata(hObject,handles);


% --- Executes on button press in pushbutton2_select_surface_motion.
function pushbutton2_select_surface_motion_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_select_surface_motion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select free surface motion(s)...';
[motion_file_name,motion_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','on');

if ~isequal(motion_dir_name,0)
    start_dir0 = motion_dir_name;
end

surface_motion = load(fullfile(motion_dir_name,motion_file_name));
plotMotion(surface_motion,'unknown',1,motion_file_name,0);

handles.metricdata.motion_file_name = motion_file_name;
handles.metricdata.motion_dir_name = motion_dir_name;

handles.metricdata.select_motion_complete = 1;

guidata(hObject,handles);



% --- Executes on button press in pushbutton3_go.
function pushbutton3_go_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.metricdata.load_complete == 0
    msgbox('Please select a profile.','Error...');
elseif handles.metricdata.select_motion_complete == 0
    msgbox('Please select a free surface motion.','Error...');
elseif handles.metricdata.load_complete == -1 % single rock layer
    h = handles.metricdata.thickness;
    Vs = handles.metricdata.Vs;
    rho = handles.metricdata.density;
    xi = handles.metricdata.damping;
    
    profile = [h,Vs,xi,rho;0,Vs,xi,rho];
    motion_file_name = handles.metricdata.motion_file_name;
    motion_dir_name = handles.metricdata.motion_dir_name;

    if ~iscell(motion_file_name) % only one motion was chosen
        fprintf('Deconvolving... ');
        
        motion = importdata(fullfile(motion_dir_name,motion_file_name));
        
        if strcmp(handles.metricdata.type_of_motion_desired,'incident')
            response = linearSiteResp(profile,motion,'elastic',0,'y'); % deconvolution
            response(:,2) = response(:,2)/2;
        elseif strcmp(handles.metricdata.type_of_motion_desired,'rock-outcrop')
            response = linearSiteResp(profile,motion,'elastic',0,'y'); % deconvolution
        elseif strcmp(handles.metricdata.type_of_motion_desired,'total')
            response = linearSiteResp(profile,motion,'rigid',0,'y'); % deconvolution
        end
        
        fh = figure;
        width = 5; height = 3;
        set(fh,'units','inches','position',[4+(5-width)/2,2+(4-height)/2,width,height]);
        % subplot(211); 
        plot(motion(:,1),motion(:,2),'b');
        hold on;
        % ylabel('Free surface motion','fontsize',12);
        %subplot(212); 
        plot(response(:,1),response(:,2),'r');
        xlabel('Time [sec]','fontsize',12);
        ylabel('Acceleration','fontsize',12);
        grid on;
        legend('Free surface motion','Deconvoled motion','location','best');
        
        [~,str1,str2] = fileparts(motion_file_name);
        new_filename = sprintf('%s_deconvolved%s',str1,str2);
        dlmwrite(fullfile(motion_dir_name,new_filename),response,'delimiter','\t','precision',6,'newline','pc');
        
        fprintf(' done.\n');
        if ispc()
            command_text = sprintf('explorer.exe %s',motion_dir_name);
            system(command_text);
        end
    else % more than one motion was chosen
        fprintf('Deconvolving... ');
        nn = length(motion_file_name);
        for sj = 1 : 1 : nn
            motion = importdata(fullfile(motion_dir_name,motion_file_name{sj}));
            
            if strcmp(handles.metricdata.type_of_motion_desired,'incident')
                response = linearSiteResp(profile,motion,'elastic',0,'y'); % deconvolution
                response(:,2) = response(:,2)/2;
            elseif strcmp(handles.metricdata.type_of_motion_desired,'rock-outcrop')
                response = linearSiteResp(profile,motion,'elastic',0,'y'); % deconvolution
            elseif strcmp(handles.metricdata.type_of_motion_desired,'total')
                response = linearSiteResp(profile,motion,'rigid',0,'y'); % deconvolution
            end
            
            fh = figure;
            width = 5; height = 3;
            set(fh,'units','inches','position',[4+(5-width)/2,2+(4-height)/2,width,height]);
            % subplot(211); 
            plot(motion(:,1),motion(:,2),'b');
            hold on;
            % ylabel('Free surface motion','fontsize',12);
            %subplot(212); 
            plot(response(:,1),response(:,2),'r');
            xlabel('Time [sec]','fontsize',12);
            ylabel('Acceleration','fontsize',12);
            grid on;
            legend('Free surface motion','Deconvoled motion','location','best');
            
            [~,str1,str2] = fileparts(motion_file_name{sj});
            new_filename = sprintf('%s_deconvolved%s',str1,str2);
            dlmwrite(fullfile(motion_dir_name,new_filename),response,'delimiter','\t','precision',6,'newline','pc');
        end
        fprintf(' done.\n');
        if ispc()
            command_text = sprintf('explorer.exe %s',motion_dir_name);
            system(command_text);
        end
    end
elseif handles.metricdata.load_complete == 1 % multiple soil layers
    profile_fname = handles.metricdata.profile_file_name;
    profile_dname = handles.metricdata.profile_dir_name;
    profile = importdata(fullfile(profile_dname,profile_fname));
    
    motion_file_name = handles.metricdata.motion_file_name;
    motion_dir_name = handles.metricdata.motion_dir_name;
    
    if ~iscell(motion_file_name) % if only one motion was chosen
        fprintf('Deconvolving... ');
        motion = importdata(fullfile(motion_dir_name,motion_file_name));
        
        if strcmp(handles.metricdata.type_of_motion_desired,'incident')
            response = linearSiteResp(profile,motion,'elastic',0,'y'); % deconvolution
            response(:,2) = response(:,2)/2;
        elseif strcmp(handles.metricdata.type_of_motion_desired,'rock-outcrop')
            response = linearSiteResp(profile,motion,'elastic',0,'y'); % deconvolution
        elseif strcmp(handles.metricdata.type_of_motion_desired,'total')
            response = linearSiteResp(profile,motion,'rigid',0,'y'); % deconvolution
        end
        
        fh = figure;
        width = 5; height = 3;
        set(fh,'units','inches','position',[4+(5-width)/2,2+(4-height)/2,width,height]);
        % subplot(211); 
        plot(motion(:,1),motion(:,2),'b');
        hold on;
        % ylabel('Free surface motion','fontsize',12);
        %subplot(212); 
        plot(response(:,1),response(:,2),'r');
        xlabel('Time [sec]','fontsize',12);
        ylabel('Acceleration','fontsize',12);
        grid on;
        legend('Free surface motion','Deconvoled motion','location','best');
        
        [~,str1,str2] = fileparts(motion_file_name);
        new_filename = sprintf('%s_deconvolved%s',str1,str2);
        dlmwrite(fullfile(motion_dir_name,new_filename),response,'delimiter','\t','precision',6,'newline','pc');
        
        fprintf(' done.\n');
        if ispc()
            command_text = sprintf('explorer.exe %s',motion_dir_name);
            system(command_text);
        end
    else % if more than one motion was chosen
        fprintf('Deconvolving... ');
        nn = length(motion_file_name);
        for sj = 1 : 1 : nn
            motion = importdata(fullfile(motion_dir_name,motion_file_name{sj}));
            
            if strcmp(handles.metricdata.type_of_motion_desired,'incident')
                response = linearSiteResp(profile,motion,'elastic',0,'y'); % deconvolution
                response(:,2) = response(:,2)/2;
            elseif strcmp(handles.metricdata.type_of_motion_desired,'rock-outcrop')
                response = linearSiteResp(profile,motion,'elastic',0,'y'); % deconvolution
            elseif strcmp(handles.metricdata.type_of_motion_desired,'total')
                response = linearSiteResp(profile,motion,'rigid',0,'y'); % deconvolution
            end
            
            fh = figure;
            width = 5; height = 3;
            set(fh,'units','inches','position',[4+(5-width)/2,2+(4-height)/2,width,height]);
            % subplot(211); 
            plot(motion(:,1),motion(:,2));
            hold on;
            % ylabel('Free surface motion','fontsize',12);
            %subplot(212); 
            plot(response(:,1),response(:,2));
            xlabel('Time [sec]','fontsize',12);
            ylabel('Acceleration','fontsize',12);
            grid on;
            legend('Free surface motion','Deconvoled motion','location','best');
            
            [~,str1,str2] = fileparts(motion_file_name{sj});
            new_filename = sprintf('%s_deconvolved%s',str1,str2);
            dlmwrite(fullfile(motion_dir_name,new_filename),response,'delimiter','\t','precision',6,'newline','pc');
        end
        
        fprintf(' done.\n');
        if ispc()
            command_text = sprintf('explorer.exe %s',motion_dir_name);
            system(command_text);
        end
    end
    
end

% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1' % uniform rock layer
        set(handles.metricdata.damping_handle,'enable','on');
        set(handles.metricdata.Vs_handle,'enable','on');
        set(handles.metricdata.density_handle,'enable','on');
        set(handles.metricdata.thickness_handle,'enable','on');
        set(handles.metricdata.load_profile_handle,'enable','off');
        handles.metricdata.load_complete = -1; % no need to load soil profile
    case 'radiobutton2' % multiple soil layers
        set(handles.metricdata.damping_handle,'enable','off');
        set(handles.metricdata.Vs_handle,'enable','off');
        set(handles.metricdata.density_handle,'enable','off');
        set(handles.metricdata.thickness_handle,'enable','off');
        set(handles.metricdata.load_profile_handle,'enable','on');
        handles.metricdata.load_complete = 0;
end
guidata(hObject,handles);


% --- Executes on button press in pushbutton4_load_profile.
function pushbutton4_load_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_load_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select profile...';
[profile_file_name,profile_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');

if ~isequal(profile_dir_name,0)
    start_dir0 = profile_dir_name;
end

Vs_profile = load(fullfile(profile_dir_name,profile_file_name));
plotVsProfileFromMatrix(Vs_profile,'y');
title(profile_file_name,'interpreter','none');

handles.metricdata.profile_file_name = profile_file_name;
handles.metricdata.profile_dir_name = profile_dir_name;

handles.metricdata.load_complete = 1;

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton4_load_profile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton4_load_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.load_profile_handle = hObject;
set(handles.metricdata.load_profile_handle,'enable','off');
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton3'
        handles.metricdata.type_of_motion_desired = 'incident';
    case 'radiobutton4'
        handles.metricdata.type_of_motion_desired = 'rock-outcrop';
    case 'radiobutton5'
        handles.metricdata.type_of_motion_desired = 'total';
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uipanel4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.metricdata.type_of_motion_desired = 'incident';
guidata(hObject,handles);


% --- Executes on button press in pushbutton5_close_all.
function pushbutton5_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;
