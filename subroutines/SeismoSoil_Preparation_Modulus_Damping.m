function varargout = SeismoSoil_Preparation_Modulus_Damping(varargin)
% SEISMOSOIL_PREPARATION_MODULUS_DAMPING MATLAB code for SeismoSoil_Preparation_Modulus_Damping.fig
%      SEISMOSOIL_PREPARATION_MODULUS_DAMPING, by itself, creates a new SEISMOSOIL_PREPARATION_MODULUS_DAMPING or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_PREPARATION_MODULUS_DAMPING returns the handle to a new SEISMOSOIL_PREPARATION_MODULUS_DAMPING or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_PREPARATION_MODULUS_DAMPING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_PREPARATION_MODULUS_DAMPING.M with the given input arguments.
%
%      SEISMOSOIL_PREPARATION_MODULUS_DAMPING('Property','Value',...) creates a new SEISMOSOIL_PREPARATION_MODULUS_DAMPING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Preparation_Modulus_Damping_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Preparation_Modulus_Damping_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Preparation_Modulus_Damping

% Last Modified by GUIDE v2.5 08-Sep-2017 01:02:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Preparation_Modulus_Damping_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Preparation_Modulus_Damping_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Preparation_Modulus_Damping is made visible.
function SeismoSoil_Preparation_Modulus_Damping_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Preparation_Modulus_Damping (see VARARGIN)

% Choose default command line output for SeismoSoil_Preparation_Modulus_Damping
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Preparation_Modulus_Damping wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Preparation_Modulus_Damping_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes during object creation, after setting all properties.
function pushbutton1_select_profile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1_select_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.select_profile_complete = 0;

guidata(hObject,handles);



% --- Executes on button press in pushbutton1_select_profile.
function pushbutton1_select_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_select_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select Vs profile data file...';
[profile_file_name,profile_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');

if ~isequal(profile_dir_name,0)
    start_dir0 = profile_dir_name;
end

soil_profile = importdata(fullfile(profile_dir_name,profile_file_name));
plotVsProfileFromMatrix(soil_profile);

checkInputs(soil_profile,'vs_profile');

handles.metricdata.profile_file_name = profile_file_name;
handles.metricdata.profile_dir_name = profile_dir_name;
handles.metricdata.soil_profile = soil_profile;

if findstr('profile_',profile_file_name) % if the profile filename contains "profile_"
    sitecode = profile_file_name(9:end-4);
else
    sitecode = profile_file_name;
end
handles.metricdata.sitecode = sitecode;
set(handles.edit1_profile_name,'string',sitecode);

handles.metricdata.select_profile_complete = 1;

guidata(hObject,handles);



function edit1_profile_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_profile_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_profile_name as text
%        str2double(get(hObject,'String')) returns contents of edit1_profile_name as a double

handles.metricdata.sitecode = get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit1_profile_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_profile_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function pushbutton2_calc_curves_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2_calc_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.calculate_curves_complete = 0;
guidata(hObject,handles);


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

handles.metricdata.PI = str2double(get(hObject,'String'));
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.metricdata.PI = 20;
guidata(hObject,handles);



% --- Executes on button press in pushbutton2_calc_curves.
function pushbutton2_calc_curves_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_calc_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;

% PI = handles.metricdata.PI;  % commented out on 3/27/2016

% strain = [0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1 3]';
% strain = [0.0001 0.0003 0.0006 0.001 0.003 0.006 ...
%         0.01 0.03 0.06 0.1 0.3 0.6 1 3 6]';  % unit: percent

strain = logspace(-4,0.778151250383644,50)'; % unit: percent

if handles.metricdata.select_profile_complete == 1

    A = handles.metricdata.soil_profile;
    sitecode = handles.metricdata.sitecode;
    profile_dir_name = handles.metricdata.profile_dir_name;
    
    nr_layer = size(A,1)-1; % number of "true" layers, i.e., excluding the "dummy" layer at the last
    
    h = A(1:end-1,1); % layer thicknesses
    Vs = A(1:end-1,2);
    rho = A(1:end-1,4); % mass density of soil of each layer
    
    PI = zeros(length(Vs),1);
    for j = 1 : 1 : length(Vs)
        if Vs(j) <= 200
            PI(j) = 10;
        elseif Vs(j) <= 360
            PI(j) = 5;
        else
            PI(j) = 0;
        end
    end
    fprintf('PI = ');
    fprintf('%d, ',PI);
    fprintf('\n');
    
    curve_mtrx = zeros(length(strain),4*nr_layer); % preallocation
    
    sigma_v0 = computeVerticalStress(h,rho); % unit of stress: Pa
    [G_over_Gmax,D,gamma_r] = darendeli2001(strain/100,sigma_v0,PI);
    
    for j = 1 : 1 : nr_layer
        curve_mtrx(:,(j-1)*4+1) = strain;
        curve_mtrx(:,(j-1)*4+3) = strain;
        curve_mtrx(:,(j-1)*4+2) = G_over_Gmax(:,j);
        curve_mtrx(:,(j-1)*4+4) = D(:,j)*100;
    end
    
    out_filename = fullfile(profile_dir_name,sprintf('curve_%s.txt',sitecode));
    dlmwrite(out_filename,curve_mtrx,'delimiter','\t','precision',6);
    showFileInExplorer(out_filename);
    
    handles.metricdata.curve_mtrx = curve_mtrx;
    handles.metricdata.calculate_curves_complete = 1;
    
    plotCurves(curve_mtrx);
       
    guidata(hObject,handles);

else
    msgbox('You need to select a Vs profile file first.','Error...');
end



% --- Executes during object creation, after setting all properties.
function checkbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.metricdata.showfig = 1;

guidata(hObject,handles);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
showfig = get(hObject,'Value');
handles.metricdata.showfig = showfig;

guidata(hObject,handles);


% --- Executes on button press in pushbutton4_H2_fitting.
function pushbutton4_H2_fitting_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_H2_fitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;

h_running = msgbox('Curve-fitting in progess. Please do not click other buttons.','Message');

if handles.metricdata.calculate_curves_complete == 1
    curve = handles.metricdata.curve_mtrx;
    H2n = gaH2(curve,handles.metricdata.showfig);
    
    sitecode = handles.metricdata.sitecode;
    profile_dir_name = handles.metricdata.profile_dir_name;
    
    out_file_name = fullfile(profile_dir_name,sprintf('H2_n_%s.txt',sitecode));
    dlmwrite(out_file_name,H2n,'delimiter','\t','precision',6);
    showFileInExplorer(out_file_name);
    
else
    msgbox('You need to calculate G/Gmax and damping curves first.','Error...');
end

if ishandle(h_running) % if the users haven't close the figure
    close(h_running);
end
msgbox('H2 curve-fitting complete!','Finished');


% --- Executes on button press in pushbutton5_H4_fitting.
function pushbutton5_H4_fitting_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_H4_fitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;

h_running = msgbox('Curve-fitting in progess. Please do not click other buttons.','Message');

if handles.metricdata.calculate_curves_complete == 1
    curve = handles.metricdata.curve_mtrx;
    H4G = fitMKZ(curve,handles.metricdata.showfig);
    H4x = gaH4x(curve,handles.metricdata.showfig);
    
    sitecode = handles.metricdata.sitecode;
    profile_dir_name = handles.metricdata.profile_dir_name;
    dlmwrite(fullfile(profile_dir_name,sprintf('H4_G_%s.txt',sitecode)),H4G,'delimiter','\t','precision',6);
    dlmwrite(fullfile(profile_dir_name,sprintf('H4_x_%s.txt',sitecode)),H4x,'delimiter','\t','precision',6);
    showFileInExplorer(fullfile(profile_dir_name,sprintf('H4_x_%s.txt',sitecode)));
else
    msgbox('You need to calculate G/Gmax and damping curves first.','Error...');
end

if ishandle(h_running) % if the users haven't close the figure
    close(h_running);
end
msgbox('H4 curve-fitting complete!','Finished');


% --- Executes on button press in pushbutton6_close_all.
function pushbutton6_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;


% --- Executes on button press in pushbutton7_return.
function pushbutton7_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Preparation_Modulus_Damping;
SeismoSoil_Preparation;
