function varargout = SeismoSoil_Preparation_Vs_Profile(varargin)
% SEISMOSOIL_PREPARATION_VS_PROFILE MATLAB code for SeismoSoil_Preparation_Vs_Profile.fig
%      SEISMOSOIL_PREPARATION_VS_PROFILE, by itself, creates a new SEISMOSOIL_PREPARATION_VS_PROFILE or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_PREPARATION_VS_PROFILE returns the handle to a new SEISMOSOIL_PREPARATION_VS_PROFILE or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_PREPARATION_VS_PROFILE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_PREPARATION_VS_PROFILE.M with the given input arguments.
%
%      SEISMOSOIL_PREPARATION_VS_PROFILE('Property','Value',...) creates a new SEISMOSOIL_PREPARATION_VS_PROFILE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Preparation_Vs_Profile_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Preparation_Vs_Profile_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Preparation_Vs_Profile

% Last Modified by GUIDE v2.5 08-Sep-2017 00:53:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Preparation_Vs_Profile_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Preparation_Vs_Profile_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Preparation_Vs_Profile is made visible.
function SeismoSoil_Preparation_Vs_Profile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Preparation_Vs_Profile (see VARARGIN)

% Choose default command line output for SeismoSoil_Preparation_Vs_Profile
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Preparation_Vs_Profile wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Preparation_Vs_Profile_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function uitable2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'data',cell(1,5));
handles.metricdata.profile_table_handle = hObject;
handles.metricdata.soil_property_complete = 0;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit1_nr_of_soil_layers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_nr_of_soil_layers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.metricdata.nr_soil_layers = 1;
guidata(hObject,handles);



function edit1_nr_of_soil_layers_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_nr_of_soil_layers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_nr_of_soil_layers as text
%        str2double(get(hObject,'String')) returns contents of edit1_nr_of_soil_layers as a double


nr_soil_layers = str2double(get(hObject,'String'));
set(handles.metricdata.profile_table_handle,'data',cell(nr_soil_layers,5));

handles.metricdata.nr_soil_layers = nr_soil_layers;
guidata(hObject,handles);


% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

profile_without_rock = get(hObject,'data');
handles.metricdata.profile_without_rock = profile_without_rock;
handles.metricdata.soil_property_complete = 1;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uipanel1_density_or_unit_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1_density_or_unit_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.density_or_unit_weight = 'density';
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel1_density_or_unit_weight.
function uipanel1_density_or_unit_weight_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1_density_or_unit_weight 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1a'
        density_or_unit_weight = 'density';
    case 'radiobutton1b'
        density_or_unit_weight = 'unit_weight';
end
handles.metricdata.density_or_unit_weight = density_or_unit_weight;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uipanel2_unit_of_damping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2_unit_of_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.unit_of_damping = 'percent';
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel2_unit_of_damping.
function uipanel2_unit_of_damping_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2_unit_of_damping 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton2a'
        unit_of_damping = 'unity';
    case 'radiobutton2b'
        unit_of_damping = 'percent';
end
handles.metricdata.unit_of_damping = unit_of_damping;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uitable3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'data',{'    Infinity',[],[],[],'   N/A'});
handles.metricdata.rock_table_handle = hObject;

handles.metricdata.rock_property_complete = 0;
guidata(hObject,handles);


% --- Executes when entered data in editable cell(s) in uitable3.
function uitable3_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable3 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

rock_property = get(hObject,'data');
handles.metricdata.rock_property = rock_property;
% disp(rock_property);
handles.metricdata.rock_property_complete = 1;
guidata(hObject,handles);


% --- Executes on button press in pushbutton1_plot_Vs_profile.
function pushbutton1_plot_Vs_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_plot_Vs_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

g = 9.81;

if handles.metricdata.rock_property_complete * handles.metricdata.soil_property_complete == 1
    density_or_unit_weight = handles.metricdata.density_or_unit_weight;
    unit_of_damping = handles.metricdata.unit_of_damping;
    
    %profile_without_rock = cell2mat(handles.metricdata.profile_without_rock);
    cell_profile = handles.metricdata.profile_without_rock;
    nrow = size(cell_profile,1);
    ncol = size(cell_profile,2);
    profile_without_rock = zeros(nrow, ncol);
    for j = 1 : nrow
        for k = 1 : ncol
            if ischar(cell_profile{j,k})
                profile_without_rock(j,k) = str2double(cell_profile{j,k});
            elseif isfloat(cell_profile{j,k})
                profile_without_rock(j,k) = cell_profile{j,k};
            else
                msgbox(sprintf('Data at row #%d column #d is invalid!',j,k),'Error');
            end
        end
    end
    
    rock_property = handles.metricdata.rock_property; % this is a 1x5 cell array
    rock_Vs = rock_property{2};
    rock_xi = rock_property{3};
    rock_rho = rock_property{4};
    if ischar(rock_Vs)
        rock_Vs = str2double(rock_Vs);
    end
    if ischar(rock_xi)
        rock_xi = str2double(rock_xi);
    end
    if ischar(rock_rho)
        rock_rho = str2double(rock_rho);
    end
    
    if isempty(rock_Vs), rock_Vs = 0; end
    if isempty(rock_xi), rock_xi = 0; end
    if isempty(rock_rho), rock_rho = 0; end
    
    if size(profile_without_rock,2) <= 1
        msgbox({'You haven''t finished entering all the properties.'},...
        'Error...');
    elseif size(profile_without_rock,2) == 2
        rr = size(profile_without_rock,1);
        profile_without_rock = [profile_without_rock,zeros(rr,3)];
        vs_profile = [profile_without_rock;0,rock_Vs,rock_xi,rock_rho,0];
        if strcmp(unit_of_damping,'percent')
            vs_profile(:,3) = vs_profile(:,3)/100;
        end
        if strcmp(density_or_unit_weight,'unit_weight')
            vs_profile(:,4) = vs_profile(:,4)*1000/g;
        end
        plotVsProfileFromMatrix(vs_profile);
    else
        vs_profile = [profile_without_rock;0,rock_Vs,rock_xi,rock_rho,0];
        if strcmp(unit_of_damping,'percent')
            vs_profile(:,3) = vs_profile(:,3)/100;
        end
        if strcmp(density_or_unit_weight,'unit_weight')
            vs_profile(:,4) = vs_profile(:,4)*1000/g;
        end
        plotVsProfileFromMatrix(vs_profile);
    end
else
    msgbox({'You haven''t finished entering all the properties.'},...
        'Error...');
end


% --- Executes on button press in pushbutton2_save.
function pushbutton2_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;

g = 9.81;

if handles.metricdata.rock_property_complete * handles.metricdata.soil_property_complete == 1
    density_or_unit_weight = handles.metricdata.density_or_unit_weight;
    unit_of_damping = handles.metricdata.unit_of_damping;
    
    cell_profile = handles.metricdata.profile_without_rock;
    nrow = size(cell_profile,1);
    ncol = size(cell_profile,2);
    profile_without_rock = zeros(nrow, ncol);
    for j = 1 : nrow
        for k = 1 : ncol
            if ischar(cell_profile{j,k})
                profile_without_rock(j,k) = str2double(cell_profile{j,k});
            elseif isfloat(cell_profile{j,k})
                profile_without_rock(j,k) = cell_profile{j,k};
            else
                msgbox(sprintf('Data at row #%d column #d is invalid!',j,k),'Error');
            end
        end
    end
    
    rock_property = handles.metricdata.rock_property; % this is a 1x5 cell array
    rock_Vs = rock_property{2};
    rock_xi = rock_property{3};
    rock_rho = rock_property{4};
    if ischar(rock_Vs)
        rock_Vs = str2double(rock_Vs);
    end
    if ischar(rock_xi)
        rock_xi = str2double(rock_xi);
    end
    if ischar(rock_rho)
        rock_rho = str2double(rock_rho);
    end
    
    vs_profile = [profile_without_rock;0,rock_Vs,rock_xi,rock_rho,0];
    
    if strcmp(unit_of_damping,'percent')
        vs_profile(:,3) = vs_profile(:,3)/100;
    end
    if strcmp(density_or_unit_weight,'unit_weight')
        vs_profile(:,4) = vs_profile(:,4)*1000/g;
    end
    
    [filename,pathname] = uiputfile('profile.txt','Save Vs profile as...');
    dlmwrite(fullfile(pathname,filename),vs_profile,'delimiter','\t');
else
    msgbox({'You haven''t finished entering all the properties.'},...
        'Error...');
end


% --- Executes on button press in pushbutton4_auto_calc_xi_rho.
function pushbutton4_auto_calc_xi_rho_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_auto_calc_xi_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

unit_of_damping = handles.metricdata.unit_of_damping;

if handles.metricdata.rock_property_complete * handles.metricdata.soil_property_complete == 0
    msgbox('You need to enter both soil and rock Vs info.','Warning');
else
    rock_property = handles.metricdata.rock_property;
    rock_Vs = rock_property{2};
    
    soil_Vs = zeros(size(handles.metricdata.profile_without_rock,1),1);
    for j = 1 : length(soil_Vs)
        soil_Vs(j) = str2double(handles.metricdata.profile_without_rock{j,2});
    end
    soil_xi = zeros(size(soil_Vs));
    soil_rho = zeros(size(soil_xi));
    nr_soil_layer = size(soil_Vs,1);

    switch unit_of_damping
        case 'unity'
            damping_unit_factor = 1;
        case 'percent'
            damping_unit_factor = 100;
    end

    if rock_Vs < 250
        rock_xi = 0.05*damping_unit_factor;
        handles.metricdata.rock_property{3} = rock_xi;
    elseif rock_Vs < 750
        rock_xi = 0.02*damping_unit_factor;
        handles.metricdata.rock_property{3} = rock_xi;
    else
        rock_xi = 0.01*damping_unit_factor;
        handles.metricdata.rock_property{3} = rock_xi;
    end

    if rock_Vs < 200
        rock_rho = 1600;
        handles.metricdata.rock_property{4} = rock_rho;
    elseif rock_Vs < 800
        rock_rho = 1800;
        handles.metricdata.rock_property{4} = rock_rho;
    else
        rock_rho = 2000;
        handles.metricdata.rock_property{4} = rock_rho;
    end

    for i = 1 : 1 : nr_soil_layer
        if soil_Vs(i) < 250
            soil_xi(i) = 0.05*damping_unit_factor;
            handles.metricdata.profile_without_rock{i,3} = soil_xi(i);
        elseif soil_Vs(i) < 750
            soil_xi(i) = 0.02*damping_unit_factor;
            handles.metricdata.profile_without_rock{i,3} = soil_xi(i);
        else
            soil_xi(i) = 0.01*damping_unit_factor;
            handles.metricdata.profile_without_rock{i,3} = soil_xi(i);
        end

        if soil_Vs(i) < 200
            soil_rho(i) = 1600;
            handles.metricdata.profile_without_rock{i,4} = soil_rho(i);
        elseif soil_Vs(i) < 800
            soil_rho(i) = 1800;
            handles.metricdata.profile_without_rock{i,4} = soil_rho(i);
        else
            soil_rho(i) = 2000;
            handles.metricdata.profile_without_rock{i,4} = soil_rho(i);
        end
        handles.metricdata.profile_without_rock{i,5} = i;
    end

    set(handles.metricdata.profile_table_handle,'data',handles.metricdata.profile_without_rock);

    set(handles.metricdata.rock_table_handle,'data',...
        {'    Infinity',handles.metricdata.rock_property{2},...
        handles.metricdata.rock_property{3},...
                handles.metricdata.rock_property{4},'   N/A'});

    guidata(hObject,handles);
end

% --- Executes on button press in pushbutton3_return.
function pushbutton3_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close SeismoSoil_Preparation_Vs_Profile;
SeismoSoil_Preparation;
