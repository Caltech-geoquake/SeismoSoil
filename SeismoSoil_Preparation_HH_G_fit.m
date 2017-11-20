function varargout = SeismoSoil_Preparation_HH_G_fit(varargin)
% SEISMOSOIL_PREPARATION_HH_G_FIT MATLAB code for SeismoSoil_Preparation_HH_G_fit.fig
%      SEISMOSOIL_PREPARATION_HH_G_FIT, by itself, creates a new SEISMOSOIL_PREPARATION_HH_G_FIT or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_PREPARATION_HH_G_FIT returns the handle to a new SEISMOSOIL_PREPARATION_HH_G_FIT or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_PREPARATION_HH_G_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_PREPARATION_HH_G_FIT.M with the given input arguments.
%
%      SEISMOSOIL_PREPARATION_HH_G_FIT('Property','Value',...) creates a new SEISMOSOIL_PREPARATION_HH_G_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_Preparation_HH_G_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_Preparation_HH_G_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_Preparation_HH_G_fit

% Last Modified by GUIDE v2.5 27-Sep-2017 15:59:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_Preparation_HH_G_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_Preparation_HH_G_fit_OutputFcn, ...
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


% --- Executes just before SeismoSoil_Preparation_HH_G_fit is made visible.
function SeismoSoil_Preparation_HH_G_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_Preparation_HH_G_fit (see VARARGIN)

% Choose default command line output for SeismoSoil_Preparation_HH_G_fit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeismoSoil_Preparation_HH_G_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_Preparation_HH_G_fit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_step1b_choose_GGmax_data.
function pushbutton_step1b_choose_GGmax_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_step1b_choose_GGmax_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select "curve" file...';
[curve_file_name,curve_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');
disp(['User selected ',fullfile(curve_dir_name,curve_file_name)]);
start_dir0 = curve_dir_name;
curve = importdata(fullfile(curve_dir_name,curve_file_name));

checkInputs(curve,'curve');

handles.metricdata.curve_file_name = curve_file_name;
handles.metricdata.curve_dir_name = curve_dir_name;
handles.metricdata.curve = curve;
handles.metricdata.step1b1 = 1;
handles.metricdata.step1b = handles.metricdata.step1b1 * handles.metricdata.step1b2;

plotCurves(curve);

% if handles.metricdata.step1a * handles.metricdata.step1b == 1
if handles.metricdata.step1b == 1  % damping parameter fitting does not require Vs profile
    set(handles.metricdata.handle_Start_damping,'enable','on');
end
if handles.metricdata.step2 * handles.metricdata.step1a == 1
    set(handles.metricdata.handle_Start,'enable','on');
end
if handles.metricdata.step2 * handles.metricdata.step1b * handles.metricdata.step1a == 1
    set(handles.metricdata.handle_Start,'enable','on');
end
guidata(hObject,handles);


% --- Executes on button press in pushbutton_step1a_choose_Vs_profile.
function pushbutton_step1a_choose_Vs_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_step1a_choose_Vs_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select Vs profile...';
[profile_file_name,profile_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');
disp(['User selected ',fullfile(profile_dir_name,profile_file_name)]);
start_dir0 = profile_dir_name;  % update start_dir0
vs_profile = importdata(fullfile(profile_dir_name,profile_file_name));

checkInputs(vs_profile,'vs_profile');

handles.metricdata.profile = vs_profile;
handles.metricdata.profile_file_name = profile_file_name;
handles.metricdata.profile_dir_name = profile_dir_name;
handles.metricdata.step1a = 1;
% if handles.metricdata.step2 * handles.metricdata.step1a == 1
%     set(handles.metricdata.handle_Start,'enable','on');
% end
% if handles.metricdata.step2 * handles.metricdata.step1b == 1
%     set(handles.metricdata.handle_Start,'enable','on');
% end
% if handles.metricdata.step1a * handles.metricdata.step1b == 1
%     set(handles.metricdata.handle_Start_damping,'enable','on');
% end
guidata(hObject,handles);


% --- Executes on button press in pushbutton_step1a_plot_Vs_profile.
function pushbutton_step1a_plot_Vs_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_step1a_plot_Vs_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vs_profile = handles.metricdata.profile;
h = vs_profile(:,1);
Vs = vs_profile(:,2);
plotVsProfileFromArrays(h,Vs,sum(h));
title(handles.metricdata.profile_file_name,'interpreter','none');


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1a'  % from Vs profile
        GGmax_data_source = 'fromVsProfile';
        set(handles.metricdata.handle_import_Vs_profile_A,'enable','on');
        set(handles.metricdata.handle_plot_Vs_profile,'enable','on');
        set(handles.metricdata.handle_import_curve,'enable','on');
        set(handles.metricdata.handle_import_curve,'enable','off');
        set(handles.metricdata.handle_import_Vs_profile_B,'enable','off');
    case 'radiobutton1b'  % from curve file
        GGmax_data_source = 'fromCurve';
        set(handles.metricdata.handle_import_Vs_profile_A,'enable','off');
        set(handles.metricdata.handle_plot_Vs_profile,'enable','off');
        set(handles.metricdata.handle_import_curve,'enable','on');
        set(handles.metricdata.handle_import_Vs_profile_B,'enable','on');
end
handles.metricdata.GGmax_data_source = GGmax_data_source;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uibuttongroup1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.GGmax_data_source = 'fromVsProfile';
guidata(hObject,handles);


% --- Executes when selected object is changed in uibuttongroup4.
function uibuttongroup4_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup4 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton2a'
        Tmax_data_source = 'fromLadd1991';
        set(handles.metricdata.handle_choose_Tau_max_file,'enable','off');
        % set(handles.metricdata.handle_Start,'enable','on');
        handles.metricdata.step2 = 1;
    case 'radiobutton2b'
        Tmax_data_source = 'imported';
        set(handles.metricdata.handle_choose_Tau_max_file,'enable','on');
        handles.metricdata.step2 = 0;
        % if handles.metricdata.step2_ == 0  % if previously Tau_max was not imported
        %     set(handles.metricdata.handle_Start,'enable','off');
        % end
end
handles.metricdata.Tmax_data_source = Tmax_data_source;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function pushbutton_step1b_choose_Vs_profile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_step1b_choose_Vs_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'enable','off');
handle_import_Vs_profile_B = hObject;  % obtain availability status
handles.metricdata.handle_import_Vs_profile_B = handle_import_Vs_profile_B;
handles.metricdata.step1b  = 0;
handles.metricdata.step1b2 = 0;
guidata(hObject,handles);


% --- Executes on button press in pushbutton_step1b_choose_Vs_profile.
function pushbutton_step1b_choose_Vs_profile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_step1b_choose_Vs_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select Vs profile...';
[profile_file_name,profile_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');
disp(['User selected ',fullfile(profile_dir_name,profile_file_name)]);
start_dir0 = profile_dir_name;  % update start_dir0
vs_profile = importdata(fullfile(profile_dir_name,profile_file_name));
plotVsProfileForGUI(vs_profile);

checkInputs(vs_profile,'vs_profile');

handles.metricdata.profile = vs_profile;
handles.metricdata.profile_file_name = profile_file_name;
handles.metricdata.profile_dir_name = profile_dir_name;
handles.metricdata.step1b2 = 1;
handles.metricdata.step1b = handles.metricdata.step1b1 * handles.metricdata.step1b2;

% if handles.metricdata.step1b == 1  % damping parameter fitting does not require Vs profile
%     set(handles.metricdata.handle_Start_damping,'enable','on');
% end
if handles.metricdata.step2 * handles.metricdata.step1a == 1
    set(handles.metricdata.handle_Start,'enable','on');
end
if handles.metricdata.step2 * handles.metricdata.step1b * handles.metricdata.step1a == 1
    set(handles.metricdata.handle_Start,'enable','on');
end

guidata(hObject,handles);


% --- Executes on button press in pushbutton_step2_choose_Tau_max_file.
function pushbutton_step2_choose_Tau_max_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_step2_choose_Tau_max_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;

filter_spec = {'*.dat;*.txt','Text files (*.dat,*.txt)';'*.*','All Files (*.*)'};
dlg_title = 'Select "Tau_max" file...';
[Tmax_file_name,Tmax_dir_name,filter_index] ...
    = uigetfile(filter_spec,dlg_title,start_dir0,'MultiSelect','off');
disp(['User selected ',fullfile(Tmax_dir_name,Tmax_file_name)]);
start_dir0 = Tmax_dir_name;
Tmax = importdata(fullfile(Tmax_dir_name,Tmax_file_name));

checkInputs(Tmax,'tau_max');

handles.metricdata.Tmax = Tmax;
handles.metricdata.step2 = 1;
handles.metricdata.step2_ = 1;
if handles.metricdata.step2 * handles.metricdata.step1a == 1
    set(handles.metricdata.handle_Start,'enable','on');
end
if handles.metricdata.step2 * handles.metricdata.step1b == 1
    set(handles.metricdata.handle_Start,'enable','on');
end
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function pushbutton_step1b_choose_GGmax_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_step1b_choose_GGmax_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'enable','off');
handle_import_curve = hObject;  % obtain availability status
handles.metricdata.handle_import_curve = handle_import_curve;
handles.metricdata.step1b  = 0;
handles.metricdata.step1b1 = 0;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton_step1a_choose_Vs_profile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_step1a_choose_Vs_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handle_import_Vs_profile_A = hObject;
handles.metricdata.handle_import_Vs_profile_A = handle_import_Vs_profile_A;
handles.metricdata.step1a = 0;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton_step1a_plot_Vs_profile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_step1a_plot_Vs_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handle_plot_Vs_profile = hObject;
handles.metricdata.handle_plot_Vs_profile = handle_plot_Vs_profile;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton_step2_choose_Tau_max_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_step2_choose_Tau_max_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'enable','off');
handle_choose_Tau_max_file = hObject;
handles.metricdata.handle_choose_Tau_max_file = handle_choose_Tau_max_file;
handles.metricdata.step2_ = 0;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pushbutton5_GGmax_para_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton5_GGmax_para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'enable','on');
handle_Start = hObject;
handles.metricdata.handle_Start = handle_Start;
guidata(hObject,handles);


% --- Executes on button press in pushbutton6_return.
function pushbutton6_return_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close SeismoSoil_Preparation_HH_G_fit;
SeismoSoil_Preparation;


% --- Executes on button press in pushbutton5_GGmax_para.
function pushbutton5_GGmax_para_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_GGmax_para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clc;

ok_to_proceed = 0;
if handles.metricdata.step2 == 0
    if handles.metricdata.step1a + handles.metricdata.step1b == 0  % both are 0
        str1 = ' 1 and 2';
    else
        str1 = ' 2';
    end
    warning_text = sprintf('You haven''t finished step(s)%s.\nYou need to finish it/them before\nstarting the analysis.',str1);
    h_msgbox = msgbox(warning_text,'Warning');
end
if handles.metricdata.step2 == 1
    if handles.metricdata.step1a + handles.metricdata.step1b > 0  % at least one is not zero
        ok_to_proceed = 1;
    else
        str1 = ' 1';
        warning_text = sprintf('You haven''t finished step(s)%s.\nYou need to finish it/them before\nstarting the analysis.',str1);
        h_msgbox = msgbox(warning_text,'Warning');
    end
end

if ok_to_proceed == 1

    if strcmpi(handles.metricdata.Tmax_data_source,'fromLadd1991')
        Tmax = [];
    end
    if strcmpi(handles.metricdata.Tmax_data_source,'imported')
        Tmax = handles.metricdata.Tmax;
    end

    if handles.metricdata.step1b == 1
        curve = handles.metricdata.curve;
    end
    vs_profile = handles.metricdata.profile;
    Vs = vs_profile(:,2);
    if vs_profile(end,1) == 0
        PI = zeros(length(Vs)-1,1);
    else
        PI = zeros(length(Vs),1);
    end
    for j = 1 : 1 : length(PI)  % PI = f(Vs), see Appendix of Shi & Asimaki (2017)
        if Vs(j) <= 200
            PI(j) = 10;
        elseif Vs(j) <= 360
            PI(j) = 5;
        else
            PI(j) = 0;
        end
    end

    show_fig = handles.metricdata.show_fig;
    save_fig = handles.metricdata.save_fig;

    %% Get HH_G parameters
    if (handles.metricdata.step1a == 1) && strcmpi(handles.metricdata.GGmax_data_source,'fromVsProfile')
        checkInputs({vs_profile,Tmax},'vs_profile_and_tau_max');
        para = hybridParaFromVsProfile(vs_profile,PI,Tmax,show_fig,save_fig,[],inf);
    elseif (handles.metricdata.step1b == 1) && strcmpi(handles.metricdata.GGmax_data_source,'fromCurve')
        checkInputs({vs_profile,Tmax},'vs_profile_and_tau_max');
        checkInputs({vs_profile,curve},'all_eql');
        [para,curves_expanded] = hybridParaFromCurves(vs_profile,curve,Tmax,show_fig,save_fig,inf);
    end

    %% Re-generate HH G/Gmax curves and damping curves (from Darendeli, 2001)
    if (handles.metricdata.step1a == 1) && strcmpi(handles.metricdata.GGmax_data_source,'fromCurve')
        length_array = size(curves_expanded,1);
        % strain_20 = logspace(-4,log10(10),length_array)'; % 20-point array from 1e-4 to 10 (unit: percent)
        nr_layer = size(para,2); % number of layers = number of columns of "para"
        curve_HH = zeros(length_array,nr_layer*4);  % pre-allocation, to put the "curve" matrix
        % h = vs_profile(1:end-1,1); % layer thicknesses
        % rho = vs_profile(1:end-1,4); % mass density of soil of each layer
        % stress = computeVerticalStress(h,rho); % unit of stress: Pa
        % [~,damping,~] = dynamicSoilParameter(strain_20/100,stress,PI); % unit of "damping": 1

        for j = 1 : 1 : nr_layer
            strain_from_curve = curve(:,(j-1)*4+1);
            Gmax = para(6,j);
            T_HH = tauHH(strain_from_curve/100,para(:,j));
            GGmax_HH = T_HH./(strain_from_curve/100)/Gmax;
            curve_HH(:,(j-1)*4+1) = strain_from_curve;
            curve_HH(:,(j-1)*4+2) = GGmax_HH;
            curve_HH(:,(j-1)*4+3) = curve(:,(j-1)*4+3);
            curve_HH(:,(j-1)*4+4) = curve(:,(j-1)*4+4);  % directly take from the original curve data
        end
    end

    if (handles.metricdata.step1a == 1) && strcmpi(handles.metricdata.GGmax_data_source,'fromVsProfile')
        length_array = 50;
        strain_20 = logspace(-4,log10(10),length_array)'; % 50-point array from 1e-4 to 10 (unit: percent)
        nr_layer = size(para,2); % number of layers = number of columns of "para"
        curve_HH = zeros(length_array,nr_layer*4);  % pre-allocation, to put the "curve" matrix
        h = vs_profile(1:end-1,1); % layer thicknesses
        if size(vs_profile,2) <= 2
            [~,rho] = getXiRho(vs_profile);
        else
            rho = vs_profile(1:end-1,4); % mass density of soil of each layer
        end
        stress = computeVerticalStress(h,rho); % unit of stress: Pa
        [~,damping,~] = dynamicSoilParameter(strain_20/100,stress,PI); % unit of "damping": 1

        for j = 1 : 1 : nr_layer
            Gmax = para(6,j);
            T_HH = tauHH(strain_20/100,para(:,j));
            GGmax_HH = T_HH./(strain_20/100)/Gmax;
            curve_HH(:,(j-1)*4+1) = strain_20;
            curve_HH(:,(j-1)*4+2) = GGmax_HH;
            curve_HH(:,(j-1)*4+3) = strain_20;
            curve_HH(:,(j-1)*4+4) = damping(:,j)*100; % unit conversion: from 1 to percent
        end
    end

    %% Save output files
    profile_file_name = handles.metricdata.profile_file_name;
    profile_dir_name = handles.metricdata.profile_dir_name;
    dir_out = profile_dir_name;

    if strcmpi(profile_file_name(1:7),'profile')
        sitecode = profile_file_name(9:end-4);
    else
        sitecode = profile_file_name(1:end-4);
    end

    HH_G_filename = sprintf('HH_G_%s.txt',sitecode);
    if (handles.metricdata.step1a == 1) && strcmpi(handles.metricdata.GGmax_data_source,'fromCurve')
        curve_HH_filename = sprintf('curve_HH_%s_expanded.txt',sitecode);
        dlmwrite(fullfile(dir_out,curve_HH_filename),curve_HH,'delimiter','\t','precision',6);
        fprintf('curve_HH file saved to directory: %s\n',dir_out);
    end
    if (handles.metricdata.step1a == 1) && strcmpi(handles.metricdata.GGmax_data_source,'fromVsProfile')
        curve_HH_filename = sprintf('curve_HH_%s.txt',sitecode);
        dlmwrite(fullfile(dir_out,curve_HH_filename),curve_HH,'delimiter','\t','precision',6);
        fprintf('curve_HH file saved to directory: %s\n',dir_out);
    end

    dlmwrite(fullfile(dir_out,HH_G_filename),para,'delimiter','\t','precision',6);
    fprintf('HH_G file saved to directory: %s\n',dir_out);

end


% --- Executes on button press in checkbox1_show_fig.
function checkbox1_show_fig_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1_show_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1_show_fig

handles.metricdata.show_fig = get(hObject,'value');
guidata(hObject,handles);


% --- Executes on button press in checkbox2_save_fig.
function checkbox2_save_fig_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2_save_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2_save_fig

handles.metricdata.save_fig = get(hObject,'value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uibuttongroup4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.Tmax_data_source = 'fromLadd1991';
handles.metricdata.step2 = 1;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function checkbox1_show_fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox1_show_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.show_fig = 0;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function checkbox2_save_fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox2_save_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.metricdata.save_fig = 1;
guidata(hObject,handles);


% --- Executes on button press in pushbutton8_close_all.
function pushbutton8_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;


