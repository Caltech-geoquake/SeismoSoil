function [ok_to_proceed,h_running] = runEquivLinearFreqDepFromGUI...
               (vs_profile,curve,nr_motion,motion,motion_name,output_dir,...
               factor_rho,factor_xi,unit_factor_accel,bedrock_type,...
               motion_type,fig_visible_option)
%
h_running = 0;

fz_axes = 12;
fz_title = 14;

%% Manually selecting the current working directory ("pwd") for Mac OS X
% % if strcmpi(computer,'maci64')
% %     current_dir = uigetdir('~','Please manually select the directory where SeismoSoil is located.');
% % end

%% Unit conversions -- soil profile
rho = vs_profile(:,4);
xi = vs_profile(:,3);
rho = rho ./ factor_rho;
xi = xi ./ factor_xi;
vs_profile(:,4) = rho;
vs_profile(:,3) = xi;

%% Bedrock type and input motion type processing
ok_to_proceed = 1;

if (strcmpi(bedrock_type,'rigid') && strcmpi(motion_type,'incident'))
    boundary = 'rigid';
    accel_division_factor = 1;
    outcrop_conversion_factor = 2;
end
if (strcmpi(bedrock_type,'rigid') && strcmpi(motion_type,'borehole'))
    boundary = 'rigid';
    accel_division_factor = 2; % input accel. will be divided by 2, becoming incident motion
    outcrop_conversion_factor = 2;
end
if (strcmpi(bedrock_type,'rigid') && strcmpi(motion_type,'outcrop'))
    boundary = 'rigid'; % because for rigid bedrock, outcrop motion = borehole motion
    accel_division_factor = 2; % input accel. will be divided by 2, becoming incident motion
    outcrop_conversion_factor = 2;
end

if (strcmpi(bedrock_type,'elastic') && strcmpi(motion_type,'incident'))
    boundary = 'elastic';
    accel_division_factor = 1; % A rule: the second number divided by the first number
    outcrop_conversion_factor = 2; % needs to give an "outcrop" motion
end
if (strcmpi(bedrock_type,'elastic') && strcmpi(motion_type,'borehole'))
    boundary = 'elastic';
    accel_division_factor = 2;
    outcrop_conversion_factor = 2;
    msgbox({'If the input motions are ''recorded at borhole'',',...
            'the bedrock should be rigid, rather than elastic.',...
            'Please modify your choices.'},'Error');
    ok_to_proceed = 0;
end
if (strcmpi(bedrock_type,'elastic') && strcmpi(motion_type,'outcrop'))
    boundary = 'elastic';
    accel_division_factor = 2;
    outcrop_conversion_factor = 2;
    choice = questdlg({'            Use With Caution',...
            'The actual motion at the soil-rock interface is',...
            'slightly different from the outcrop motion, so it is',...
            'recommended that the users deconvolve the rock',...
            'outcrop motion to incident motion.',...
            'Proceed anyway?'}, ...
	'Caution', ...
	'Yes','No','Yes');
    switch choice
        case 'Yes'
            ok_to_proceed = 1;
        case 'No'
            ok_to_proceed = 0;
    end
end

%% Start!
if ok_to_proceed == 1
    h_running = msgbox('Analysis in progess. Please do not click other buttons.','Calculating');
    
    mkdir(output_dir);
    
    parfor ii = 1 : nr_motion
        fprintf('Ground motion No. %d\n',ii);
        
        current_motion = motion{ii};
        [a_temp,motion_name_without_ext,ext] = fileparts(motion_name{ii});
        output_dir2 = fullfile(output_dir,motion_name_without_ext);
        mkdir(output_dir2);
        
        accel_temp = current_motion(:,2);
        accel_temp = accel_temp ./ unit_factor_accel;
        current_motion(:,2) = accel_temp;
        accel_incident = current_motion;
        accel_incident(:,2) = accel_incident(:,2)/accel_division_factor;
        
        [freq_array,tf,t_out,accel_on_surface] ...
            = callFDEQ_internal(vs_profile,accel_incident,curve,boundary,output_dir2);
        
        filename_TF = sprintf('%s_equivalent_linear_(FD)_TF%s',motion_name_without_ext,ext);
        filename_surface_accel = sprintf('%s_accel_on_surface%s',motion_name_without_ext,ext);
        % filename_new_profile = sprintf('%s_re-discretized_profile%s',motion_name_without_ext,ext);
        % filename_out_a = sprintf('%s_time_history_accel%s',motion_name_without_ext,ext);
        % filename_out_v = sprintf('%s_time_history_veloc%s',motion_name_without_ext,ext);
        % filename_out_d = sprintf('%s_time_history_displ%s',motion_name_without_ext,ext);
        % filename_out_gamma = sprintf('%s_time_history_strain%s',motion_name_without_ext,ext);
        % filename_out_tau = sprintf('%s_time_history_stress%s',motion_name_without_ext,ext);
        % filename_max_avd = sprintf('%s_max_a_v_d%s',motion_name_without_ext,ext);
        % filename_max_gt = sprintf('%s_max_gamma_tau%s',motion_name_without_ext,ext);
        
        dlmwrite(fullfile(output_dir2,filename_TF),[freq_array,tf],'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_surface_accel),[t_out,accel_on_surface],'delimiter','\t','precision',7,'newline','pc');
        % dlmwrite(fullfile(output_dir2,filename_new_profile),new_profile,'delimiter','\t','precision',6,'newline','pc');
        % dlmwrite(fullfile(output_dir2,filename_out_a),out_a,'delimiter','\t','precision',6,'newline','pc');
        % dlmwrite(fullfile(output_dir2,filename_out_v),out_v,'delimiter','\t','precision',6,'newline','pc');
        % dlmwrite(fullfile(output_dir2,filename_out_d),out_d,'delimiter','\t','precision',6,'newline','pc');
        % dlmwrite(fullfile(output_dir2,filename_out_gamma),out_gamma,'delimiter','\t','precision',6,'newline','pc');
        % dlmwrite(fullfile(output_dir2,filename_out_tau),out_tau,'delimiter','\t','precision',6,'newline','pc');
        % dlmwrite(fullfile(output_dir2,filename_max_avd),max_avd,'delimiter','\t','precision',6,'newline','pc');
        % dlmwrite(fullfile(output_dir2,filename_max_gt),max_gt,'delimiter','\t','precision',6,'newline','pc');
        
        %% Plot transfer function
        fig1 = figure('visible',fig_visible_option); % Create a figure object named "fig1" ("handle" of the fig)
        width = 5; height = 7;
        xLeft = (8-width)/2; yBottom = (10-height)/2;
        set(fig1,'Units','inches','paperposition',[xLeft,yBottom,width,height]);
        
        subplot(211);
        semilogx(freq_array,tf,'k','linewidth',1.5);
        xlabel('Frequency (Hz)','fontsize',fz_axes);
        ylabel('Amplification Factor','fontsize',fz_axes);
        grid on;
        set(gca,'fontsize',fz_axes);
        xlim([min(freq_array) max(freq_array)]);
        tf_title_text = sprintf('Frequency Dependent Equivalent Linear Amplification Factor\n%s',motion_name_without_ext);
        title(tf_title_text,'interpreter','none','fontsize',fz_title);
        
        subplot(212);
        plot(freq_array,tf,'k','linewidth',1.5);
        xlabel('Frequency (Hz)','fontsize',fz_axes);
        ylabel('Amplification Factor','fontsize',fz_axes);
        grid on;
        set(gca,'fontsize',fz_axes);
        xlim([min(freq_array) max(freq_array)]);
        
        tf_fig_filename = sprintf('%s_equivalent_linear_(FD)_TF.png',motion_name_without_ext);
        saveas(fig1,fullfile(output_dir2,tf_fig_filename));
        
        %% Plot accelerations
        fig2 = figure('visible',fig_visible_option); % Create a figure object named "fig1" ("handle" of the fig)
        width = 5; height = 4;
        xLeft = (18.5-width)/2;  yBottom = (11-height)/2;
        set(fig2,'units','inches','paperposition',[xLeft,yBottom,width,height]);
        
        plot(t_out,accel_on_surface,'color',[.6 .6 .6]); hold on;
        plot(accel_incident(:,1),accel_incident(:,2)*outcrop_conversion_factor,'r');
        ylabel('Acceleration (m/s^2)','fontsize',fz_axes);
        xlabel('Time (s)','fontsize',fz_axes);
        set(gca,'fontsize',fz_axes);
        xlim([0 max(t_out)]);
        grid on;
        h_legend = legend('Output','Rock outcrop','location','northeast');%,'fontsize',fz_axes);
        set(h_legend,'fontsize',fz_axes); % actually the fontsize of legend is by default the same as axes
        accel_title_text = sprintf('Input and Output Accelerations\n%s',motion_name_without_ext);
        title(accel_title_text,'interpreter','none','fontsize',fz_title);
        
        accel_fig_filename = sprintf('%s_accelerations.png',motion_name_without_ext);
        saveas(fig2,fullfile(output_dir2,accel_fig_filename));
        
        %% Plot max_a_v_d_gamma_tau
        % fz_axes = fz_axes - 1;
        % fz_title = fz_title - 1;
        % 
        % layer_boundary_depth = max_avd(:,1); % nr_layer x 1
        % max_a = max_avd(:,2);
        % max_v = max_avd(:,3);
        % max_d = max_avd(:,4);
        % layer_midpoint_depth = max_gt(:,1);  % (nr_layer-1) x 1
        % max_gamma = max_gt(:,2);
        % max_tau = max_gt(:,3);
        % 
        % fig3 = figure('visible',fig_visible_option);
        % width = 12; height = 7;
        % xLeft = (12.5-width)/2;  yBottom = (11-height)/2;
        % set(fig3,'units','inches','position',[xLeft,yBottom,width,height]);
        % 
        % subplot(151);
        % plot(max_a,layer_boundary_depth,'b.-');
        % set(gca,'yDir','reverse','fontsize',fz_axes);
        % ylim([0 max(layer_boundary_depth)]);
        % xlabel('Max. accel. (m/s^2)','fontsize',fz_axes);
        % ylabel('Depth (m)','fontsize',fz_axes);
        % 
        % subplot(152);
        % plot(max_v*100,layer_boundary_depth,'b.-');
        % set(gca,'yDir','reverse','fontsize',fz_axes);
        % ylim([0 max(layer_boundary_depth)]);
        % xlabel('Max. veloc. (cm/s)');
        % 
        % subplot(153);
        % plot(max_d*100,layer_boundary_depth,'b.-');
        % set(gca,'yDir','reverse','fontsize',fz_axes);
        % ylim([0 max(layer_boundary_depth)]);
        % xlabel('Max. displ. (cm)','fontsize',fz_axes);
        % title(motion_name_without_ext,'interpreter','none','fontsize',fz_title);
        % 
        % subplot(154);
        % plot(max_gamma*100,layer_midpoint_depth,'b.-');
        % set(gca,'yDir','reverse','fontsize',fz_axes);
        % ylim([0 max(layer_boundary_depth)]);
        % xlabel('Max. strain (%)','fontsize',fz_axes);
        % 
        % subplot(155);
        % plot(max_tau/1000,layer_midpoint_depth,'b.-');
        % set(gca,'yDir','reverse','fontsize',fz_axes);
        % ylim([0 max(layer_boundary_depth)]);
        % xlabel('Max. stress (kPa)','fontsize',fz_axes);
        % 
        % max_fig_filename = sprintf('%s_max_a_v_d_gamma_tau.png',motion_name_without_ext);
        % saveas(fig3,fullfile(output_dir2,max_fig_filename));
        
        %%
        if strcmpi(fig_visible_option,'off')
            close(fig1);
            close(fig2);
            % close(fig3);
        end
    end
end

end

function [freq_array,tf,t_out,accel_on_surface]...
         = callFDEQ_internal(vs_profile,motion,curve,boundary,dir_FDEQ)
%
%
% [Inputs]
%    vs_profile
%    motion
%    curve
%    fig_visible: 'off' (default) and 'on'
%    boundary: 'elastic' (default) and 'rigid'
%    output_or_not: 'n' (default) and 'y'
%
% [Outputs]
%    freq_array
%    tf: linear transfer function
%    t_out: time array, same as input time array
%    accel_on_surface: only one column
%    new_profile: new (denser) soil vs profile after re-discretization
%    out_a: acceleration time history of every layer
%    out_v
%    out_d
%    out_gamma
%    out_tau
%    max_avd: maximum acceleration, velocity, and displacement at every
%             interface (including the soil-rock interface)
%    max_gt: maximum strain (gamma) and stress (tau) of every layer
%
% NOTES:
% (1) The "motion" is incident motion (already converted to 
%     incident motion by the function that calls this function), regardless
%     of the "boundary" type.
% (2) Damping curves are adjusted, based on the initial damping values
%     provided in the "vs_profile" matrix.
%
% Jian Shi, 12/18/2013

if strcmp(boundary,'rigid')
    motion(:,2) = motion(:,2)*2; % incident --> total motion at borehole
    input_motion_type = 'borehole';
end
if strcmp(boundary,'elastic')
    motion(:,2) = motion(:,2)*2; % incident --> outcrop motion
    input_motion_type = 'outcrop';
end

current_dir = pwd;
        
[fortran_dir,~,~] = fileparts(mfilename('fullpath'));  % assume same as this M file
if ispc()
    if isdeployed()
        copyfile('FDEQ.exe',dir_FDEQ);
    else
        copyfile(fullfile(fortran_dir,'FDEQ.exe'),dir_FDEQ);
    end
elseif ismac()
    if isdeployed()
        copyfile('/Applications/SeismoSoil.app/Contents/MacOS/FDEQ',dir_FDEQ);
    else
        copyfile('FDEQ',dir_FDEQ);
    end
else  % Linux, not mac
    warning('Compiled SeismoSoil does not work properly on Linux.')
    copyfile(fullfile(fortran_dir,'FDEQ'),dir_FDEQ);
end

% If motion length > 16384, increase sampling time interval so that the
% Fortran code can perform FFT to it.
time = motion(:,1);
accel = motion(:,2);
while length(time) > 16384
    time = time(1:2:end);
    accel = accel(1:2:end);
end
motion = [time,accel];

% Adjust the damping curves according to damping values in "profile"
xi = vs_profile(:,3);
nr_layer = length(xi); % layer count, including the rock halfspace
for k = 1 : 1 : nr_layer-1
    curve(:,(k-1)*4+4) = curve(:,(k-1)*4+4) - curve(1,(k-1)*4+4) + xi(k)*100;
end

% call resp_FD, which calls FDEQ.exe
[gm,tf2col] = resp_FD_internal(input_motion_type,motion,vs_profile,curve,dir_FDEQ);

% delete the text files generated during "resp_FD"
delete('kausel001.sac');
delete('ktf001.sac');
delete('profile.dat');
delete('soil.dat');
delete('input_motion.dat');

if strcmpi(computer,'pcwin64')
    delete('FDEQ.exe');
elseif strcmpi(computer,'maci64')
    delete('./FDEQ');
end

cd(current_dir); % return to the directory before calling "resp_FD"

freq_array = tf2col(:,1);
tf = tf2col(:,2);
if freq_array(1) == 0 % delete the first element if freq(1) == 0
    freq_array = freq_array(2:end);
    tf = tf(2:end);
end
t_out = gm(:,1);
accel_on_surface = gm(:,2);

end

% % 

function [GM,transfer] = resp_FD_internal(input_accel_type,accel,vs_profile,curve,dirFD)

% function calling the FD model
%
% function [ GM, transfer] = resp_FD(boundary, acc, prof, curve, dir_FD)
%
%       GM:       surface ground motion time history
%       boundary: boundary condition
%       acc:      input motion
%       prof:     profile data
%       curve:    modulus reduction and material damping curve
%
% Edited by Jian Shi on 12/16/2013

cd(dirFD);

h = vs_profile(:,1); % height (m)
vs = vs_profile(:,2); % shear velocity (m/s)
rho = vs_profile(:,4); % density (kg/m^3)
% x=parameter(:,3); % initial value of damping ratio
mtrl_nr = vs_profile(:,5); %material number
nlayer = length(h); % include layer of bedrock

str1 = 'input_motion.dat';
str2 = 'convolution';
str3 = input_accel_type;
dlmwrite(str1,accel);

delete('profile.dat')

fid = fopen('profile.dat','a');

fprintf(fid, '%s\n', str1);
fprintf(fid, '%s\n', str2);
fprintf(fid, '%s\n', str3);
fprintf(fid, '%10.0f\n', 1);
fprintf(fid, '%10.0f %10.0f\n', 1, nlayer);
fprintf(fid, '%10.1f\n', 6.0);
fprintf(fid, '%10.1f\n', 0.1);
fprintf(fid, '%10.0f\n', 10);
fprintf(fid, '%10.0f\n', nlayer);

for k=1:nlayer
    fprintf(fid,'%10.0f %10.2f %10.3f %10.1f %10.1f\n',mtrl_nr(k),vs(k),h(k),rho(k),25.0);
end
fclose(fid);

%degradation=load('curve.dat');
degradation = curve;
n_obs = size(degradation,1);
n_ma = size(degradation,2)/4;
strn_G = zeros(n_obs, n_ma);
G_vec = zeros(n_obs, n_ma);
strn_x = zeros(n_obs, n_ma);
x_vec = zeros(n_obs, n_ma);

for k = 1:n_ma
    strn_G(:,k)=degradation(:,(k-1)*4+1);
    G_vec(:,k)=degradation(:,(k-1)*4+2);  % allocate G/Gmax according to material type
    strn_x(:,k)=degradation(:,(k-1)*4+3);
    x_vec(:,k)=degradation(:,(k-1)*4+4);  % allocate damping according to material type
end

N_obs=length(strn_G(:,1));

delete('soil.dat')
fid2 = fopen('soil.dat', 'a');

fprintf(fid2, '%6.0f\n', n_ma);

for k=1:n_ma
    fprintf(fid2, '%6.0f\n', k);
    fprintf(fid2, '%6.0f\n', N_obs);
    for j=1:N_obs
        fprintf(fid2,  '%10.4f %10.2f\n', strn_G(j,k), G_vec(j,k));
    end
    fprintf(fid2, '%6.0f\n', N_obs);
    for j=1:N_obs
        fprintf(fid2, '%10.4f %10.2f\n', strn_x(j,k), x_vec(j,k));
    end
    
end
fclose(fid2);

if strcmpi(computer,'pcwin64')
    system('FDEQ');
elseif strcmpi(computer,'maci64')
    system ('./FDEQ');
end

GM = load('kausel001.sac');
transfer = load('ktf001.sac');

end
     
     
