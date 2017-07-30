function varargout = SeismoSoil_PostProcessing(varargin)
% SEISMOSOIL_POSTPROCESSING MATLAB code for SeismoSoil_PostProcessing.fig
%      SEISMOSOIL_POSTPROCESSING, by itself, creates a new SEISMOSOIL_POSTPROCESSING or raises the existing
%      singleton*.
%
%      H = SEISMOSOIL_POSTPROCESSING returns the handle to a new SEISMOSOIL_POSTPROCESSING or the handle to
%      the existing singleton*.
%
%      SEISMOSOIL_POSTPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEISMOSOIL_POSTPROCESSING.M with the given input arguments.
%
%      SEISMOSOIL_POSTPROCESSING('Property','Value',...) creates a new SEISMOSOIL_POSTPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeismoSoil_PostProcessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeismoSoil_PostProcessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeismoSoil_PostProcessing

% Last Modified by GUIDE v2.5 09-Feb-2015 01:23:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeismoSoil_PostProcessing_OpeningFcn, ...
                   'gui_OutputFcn',  @SeismoSoil_PostProcessing_OutputFcn, ...
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


% --- Executes just before SeismoSoil_PostProcessing is made visible.
function SeismoSoil_PostProcessing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeismoSoil_PostProcessing (see VARARGIN)

% Choose default command line output for SeismoSoil_PostProcessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% When this property is set to 1, this GUI will stays open even if "close
% all" command is executed.
setappdata(hObject, 'IgnoreCloseAll', 1);

% UIWAIT makes SeismoSoil_PostProcessing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeismoSoil_PostProcessing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_select_dir.
function pushbutton1_select_dir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_select_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global start_dir0;
results_dir = uigetdir(start_dir0,'Select folder where results are stored...');
handles.metricdata.results_dir = results_dir;

file_list = listFile(results_dir,'*_accelerations.png'); % return a cell array "file_list"
accel_time_history_figure_filename = file_list{1}; % file_list should have only 1 element
result_event_name = accel_time_history_figure_filename(1:end-18);
handles.metricdata.result_event_name = result_event_name;
guidata(hObject,handles);


% --- Executes on button press in pushbutton2_plot_loops.
function pushbutton2_plot_loops_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2_plot_loops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results_dir = handles.metricdata.results_dir;
event_name = handles.metricdata.result_event_name;
strain_time_history_filename1 = fullfile(results_dir,sprintf...
    ('%s_time_history_strain.txt',event_name));
strain_time_history_filename2 = fullfile(results_dir,sprintf...
    ('%s_time_history_strain.dat',event_name));
stress_time_history_filename1 = fullfile(results_dir,sprintf...
    ('%s_time_history_stress.txt',event_name));
stress_time_history_filename2 = fullfile(results_dir,sprintf...
    ('%s_time_history_stress.dat',event_name));
profile_filename1 = fullfile(results_dir,sprintf...
    ('%s_re-discretized_profile.txt',event_name));
profile_filename2 = fullfile(results_dir,sprintf...
    ('%s_re-discretized_profile.dat',event_name));
try
    strain = importdata(strain_time_history_filename1);
    stress = importdata(stress_time_history_filename1);
    vs_profile = importdata(profile_filename1);
catch
    strain = importdata(strain_time_history_filename2);
    stress = importdata(stress_time_history_filename2);
    vs_profile = importdata(profile_filename2);
end
nr_layer = size(strain,2);
thickness = vs_profile(:,1);
depth = convertThicknessToDepth(thickness);

for j = 1 : 1 : nr_layer
    figure;
    plot(strain(:,j)*100,stress(:,j)/1000);
    xlabel('Strain [%]');
    ylabel('Stress [kPa]');
    grid on;
    title(sprintf('Layer #%d. Depth = %.2f m',j,depth(j)));
end



% --- Executes on button press in pushbutton3_plot_and_save_loops.
function pushbutton3_plot_and_save_loops_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_plot_and_save_loops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results_dir = handles.metricdata.results_dir;
event_name = handles.metricdata.result_event_name;
strain_time_history_filename1 = fullfile(results_dir,sprintf...
    ('%s_time_history_strain.txt',event_name));
strain_time_history_filename2 = fullfile(results_dir,sprintf...
    ('%s_time_history_strain.dat',event_name));
stress_time_history_filename1 = fullfile(results_dir,sprintf...
    ('%s_time_history_stress.txt',event_name));
stress_time_history_filename2 = fullfile(results_dir,sprintf...
    ('%s_time_history_stress.dat',event_name));
profile_filename1 = fullfile(results_dir,sprintf...
    ('%s_re-discretized_profile.txt',event_name));
profile_filename2 = fullfile(results_dir,sprintf...
    ('%s_re-discretized_profile.dat',event_name));
try
    strain = importdata(strain_time_history_filename1);
    stress = importdata(stress_time_history_filename1);
    vs_profile = importdata(profile_filename1);
catch
    strain = importdata(strain_time_history_filename2);
    stress = importdata(stress_time_history_filename2);
    vs_profile = importdata(profile_filename2);
end
nr_layer = size(strain,2);
thickness = vs_profile(:,1);
depth = convertThicknessToDepth(thickness);

mkdir(fullfile(results_dir,'Stress strain loops'));
for j = 1 : 1 : nr_layer
    fig = figure;
    plot(strain(:,j)*100,stress(:,j)/1000);
    grid on;
    xlabel('Strain [%]');
    ylabel('Stress [kPa]');
    title(sprintf('Layer #%d. Depth = %.2f m',j,depth(j)));
    saveas(fig,fullfile(results_dir,'Stress strain loops',...
        sprintf('Stress_strain_loops_Layer_%d.png',j)));
end

% --- Executes on button press in pushbutton4_view_movie.
function pushbutton4_view_movie_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_view_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results_dir = handles.metricdata.results_dir;
event_name = handles.metricdata.result_event_name;
displ_time_history_filename1 = fullfile(results_dir,sprintf...
    ('%s_time_history_displ.txt',event_name));
displ_time_history_filename2 = fullfile(results_dir,sprintf...
    ('%s_time_history_displ.dat',event_name));
surface_accel_filename1 = fullfile(results_dir,sprintf...
    ('%s_accel_on_surface.txt',event_name));
surface_accel_filename2 = fullfile(results_dir,sprintf...
    ('%s_accel_on_surface.dat',event_name));
profile_filename1 = fullfile(results_dir,sprintf...
    ('%s_re-discretized_profile.txt',event_name));
profile_filename2 = fullfile(results_dir,sprintf...
    ('%s_re-discretized_profile.dat',event_name));
try
    displ = importdata(displ_time_history_filename1);
    surf_accel = importdata(surface_accel_filename1);
    vs_profile = importdata(profile_filename1);
catch
    displ = importdata(displ_time_history_filename2);
    surf_accel = importdata(surface_accel_filename2);
    vs_profile = importdata(profile_filename2);
end

time = surf_accel(:,1);
dt = time(2) - time(1);
h = vs_profile(:,1);
z = convertThicknessToDepth(h);

max_displ = max(max(abs(displ)));

movie_speed = handles.metricdata.movie_speed;

figure;
for k = 1 : movie_speed : length(time)
    plot(displ(k,:)*100, z(1:end), '-.o');
    grid on;
    set(gca,'Ydir','reverse','fontsize',12);
    xlabel('Horizontal ground displacement [cm]');
    ylabel('Depth [m]');
    xlim([-max_displ*100, max_displ*100]);
    ylim([0 max(z)]);
    title(sprintf('Time = %.2f sec (%d times faster)',k*dt,movie_speed));
    pause(dt);
end


% --- Executes on button press in pushbutton5_view_and_save_movie.
function pushbutton5_view_and_save_movie_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5_view_and_save_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results_dir = handles.metricdata.results_dir;
event_name = handles.metricdata.result_event_name;
displ_time_history_filename1 = fullfile(results_dir,sprintf...
    ('%s_time_history_displ.txt',event_name));
displ_time_history_filename2 = fullfile(results_dir,sprintf...
    ('%s_time_history_displ.dat',event_name));
surface_accel_filename1 = fullfile(results_dir,sprintf...
    ('%s_accel_on_surface.txt',event_name));
surface_accel_filename2 = fullfile(results_dir,sprintf...
    ('%s_accel_on_surface.dat',event_name));
profile_filename1 = fullfile(results_dir,sprintf...
    ('%s_re-discretized_profile.txt',event_name));
profile_filename2 = fullfile(results_dir,sprintf...
    ('%s_re-discretized_profile.dat',event_name));
try
    displ = importdata(displ_time_history_filename1);
    surf_accel = importdata(surface_accel_filename1);
    vs_profile = importdata(profile_filename1);
catch
    displ = importdata(displ_time_history_filename2);
    surf_accel = importdata(surface_accel_filename2);
    vs_profile = importdata(profile_filename2);
end

time = surf_accel(:,1);
dt = time(2) - time(1);
h = vs_profile(:,1);
z = convertThicknessToDepth(h);

max_displ = max(max(abs(displ)));

movie_speed = handles.metricdata.movie_speed;  % an integer, e.g., 20

frames(floor(length(time)/movie_speed)) = struct('cdata',[],'colormap',[]);
counter = 1;
hfig = figure;
for k = 1 : movie_speed : length(time)
    plot(displ(k,:)*100, z(1:end), '-.o');
    grid on;
    set(gca,'Ydir','reverse','fontsize',12);
    xlabel('Horizontal ground displacement [cm]');
    ylabel('Depth [m]');
    xlim([-max_displ*100, max_displ*100]);
    ylim([0 max(z)]);
    title(sprintf('Time = %.2f sec (%d times faster)',k*dt,movie_speed));
    pause(dt);
    frames(counter) = getframe(hfig);
    counter = counter + 1;
end

movie_filename = fullfile(results_dir,sprintf('%s.avi',event_name));
my_video = VideoWriter(movie_filename);
my_video.FrameRate = floor(1/dt); 
% my_video.FrameRate = 50;
open(my_video);
writeVideo(my_video,frames);
close(my_video);




function edit1_relative_speed_Callback(hObject, eventdata, handles)
% hObject    handle to edit1_relative_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1_relative_speed as text
%        str2double(get(hObject,'String')) returns contents of edit1_relative_speed as a double

movie_speed = str2double(get(hObject,'String'));
if movie_speed <= 0
    msgbox('You can''t choose 0 or negative numbers.','Error');
end
handles.metricdata.movie_speed = floor(movie_speed);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit1_relative_speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1_relative_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.metricdata.movie_speed = 20;
guidata(hObject,handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
