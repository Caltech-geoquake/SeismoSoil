function [ok_to_proceed,h_running] = runNLH2FromGUI(vs_profile,curve,H2n,...
    nr_motion,motion,motion_name,output_dir,...
    factor_rho,factor_xi,unit_factor_accel,bedrock_type,motion_type,...
    fig_visible,use_fortran,use_parallel)

h_running = 0;

rho = vs_profile(:,4);
xi = vs_profile(:,3);
rho = rho ./ factor_rho;
xi = xi ./ factor_xi;
vs_profile(:,4) = rho;
vs_profile(:,3) = xi;

para = H2n;

if use_parallel == 1
    nr_cores = inf; % using a maximum of "nr_cores" workers or threads
elseif use_parallel == 0
    nr_cores = 0; % do not use multiple cores/threads
end

%% Manually selecting the current working directory ("pwd") for Mac OS X
% if strcmpi(computer,'maci64')
%     current_dir = uigetdir('~','Please manually select the directory where SeismoSoil is located.');
% end

%% Bedrock type and input motion type processing
ok_to_proceed = 1;

if (strcmpi(bedrock_type,'rigid') && strcmpi(motion_type,'incident'))
    boundary = 'rigid';
    accel_division_factor = 1; % the raw input accel. is divided by this number to participate in NL simulation
    outcrop_conversion_factor = 2; % the accel. above is multiplied by this number to make it "outcrop"
    is_incident_multiplication_factor = 2; % if "accel_in" is incident motion, this value should be 2, otherwise, 1
end
if (strcmpi(bedrock_type,'rigid') && strcmpi(motion_type,'borehole'))
    boundary = 'borehole';
    accel_division_factor = 1;
    outcrop_conversion_factor = 1;
    is_incident_multiplication_factor = 1; % if "accel_in" is incident motion, this value should be 2, otherwise, 1
end
if (strcmpi(bedrock_type,'rigid') && strcmpi(motion_type,'outcrop'))
    boundary = 'borehole'; % because for rigid bedrock, outcrop motion = borehole motion
    accel_division_factor = 1;
    outcrop_conversion_factor = 1;
    is_incident_multiplication_factor = 1; % if "accel_in" is incident motion, this value should be 2, otherwise, 1
end

if (strcmpi(bedrock_type,'elastic') && strcmpi(motion_type,'incident'))
    boundary = 'elastic';
    accel_division_factor = 1; % A rule: the second number divided by the first number
    outcrop_conversion_factor = 2; % needs to give an "outcrop" motion
    is_incident_multiplication_factor = 2; % if "accel_in" is incident motion, this value should be 2, otherwise, 1
end
if (strcmpi(bedrock_type,'elastic') && strcmpi(motion_type,'borehole'))
    boundary = 'elastic';
    accel_division_factor = 2; % input accel. will be divided by 2, becoming incident motion
    outcrop_conversion_factor = 2;
    is_incident_multiplication_factor = 2; % if "accel_in" is incident motion, this value should be 2, otherwise, 1
    msgbox({'If the input motions are ''recorded at borhole'',',...
            'the bedrock should be rigid, rather than elastic.',...
            'Please modify your choices.'},'Error');
    ok_to_proceed = 0;
end
if (strcmpi(bedrock_type,'elastic') && strcmpi(motion_type,'outcrop'))
    boundary = 'elastic';
    accel_division_factor = 2; % input accel. will be divided by 2, becoming incident motion
    outcrop_conversion_factor = 2;
    is_incident_multiplication_factor = 2; % if "accel_in" is incident motion, this value should be 2, otherwise, 1
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

if ok_to_proceed == 1
    h_running = msgbox('Analysis in progess. Please do not click other buttons.','Calculating');
    
    mkdir(output_dir);
    
    %% Public parameters
    f_max = 30; % maximum frequency modeled, unit is Hz
    ppw = 10; % points per wavelength
    n_dt = 30; % number of sub-steps in one time step
    N_spr = 120; % number of Iwan springs
    % fig_visible = 'off';
    
    fz = 14;   % fontsize
    lw = 1.75; % linewidth
    
    tabk = [...
        1.72333E-03 1.66958E-02 8.98758E-02;...
        1.80701E-03 3.81644E-02 6.84635E-02;...
        5.38887E-03 9.84666E-03 9.67052E-02;...
        1.99322E-02 -1.36803E-02 1.20172E-01;...
        8.49833E-02 -2.85125E-02 1.30728E-01;...
        4.09335E-01 -5.37309E-02 1.38746E-01;...
        2.05951E+00 -6.65035E-02 1.40705E-01;...
        1.32629E+01 -1.33696E-01 2.14647E-01];
    
    %% Public file paths
    main_directory = output_dir;
    
    %% Re-discretize soil profile
    h = vs_profile(:,1); % thickness of each layer (m)
    vs = vs_profile(:,2); % shear wave velocity of each layer (m/s)
    D = vs_profile(:,3); % initial value of damping ratio (unit: 1)
    rho = vs_profile(:,4); % mass density of soil of each layer
    material_nr = vs_profile(:,5); % material number
    [h,vs,D,rho,material_nr] = stratify(h,vs,D,rho,material_nr);
    
    new_profile = [h,vs,D,rho,material_nr];
    
    parfor (ii = 1 : nr_motion, nr_cores)
        fprintf('Ground motion No. %d\n',ii);
        %% Prepare output folder
        current_motion_name_with_ext = motion_name{ii};
        [temp_dir,motion_name_without_ext,ext] = fileparts(current_motion_name_with_ext);
        output_dir2 = fullfile(output_dir,motion_name_without_ext);
        dir_h2rev = output_dir2;
        mkdir(dir_h2rev);
        if strcmpi(computer,'pcwin64')
            if use_fortran == 1
                copyfile('NLH2.exe',dir_h2rev);
            end
        elseif strcmpi(computer,'maci64')
            % fortran_exec_dir = fullfile(current_dir,'SeismoSoil.app/Contents/MacOS/');
            if use_fortran == 1
                % copyfile(fullfile(fortran_exec_dir,'NLH2'),dir_h2rev);
                if isdeployed  % running as compiled executable
                    copyfile('/Applications/SeismoSoil.app/Contents/MacOS/NLH2',dir_h2rev);
                else  % running in MATLAB
                    copyfile('NLH2',dir_h2rev);
                end
            end
        end
        
        %% Reconstruct "tabk.dat"
        dlmwrite(fullfile(dir_h2rev,'tabk.dat'),tabk,'delimiter','\t');
        
        %% Prepare input motions
        current_motion = motion{ii};
        accel_temp = current_motion(:,2);
        accel_temp = accel_temp ./ unit_factor_accel;
        current_motion(:,2) = accel_temp;
        accel_in = current_motion;
        accel_in(:,2) = accel_in(:,2)/accel_division_factor;
        
        %% Compute central frequency of incident motion
        %     cf = centralFreq(accel_in);
        %     dlmwrite(fullfile(output_dir2,'central_freq.dat'),cf);
        
        %% Determine n_dt from the input motion
% %         dt_cr = 0.5/f_max/ppw; % the factor can range from 0.5 to 0.9
% %         dt = accel_in(2,1)-accel_in(1,1);
% %         n_dt = ceil(dt/dt_cr)*2; % times 2 for extra safety
        
        %% Doing non-linear analysis
% %         if strcmpi(computer,'pcwin64')
            current_dir = pwd;
% %         elseif strcmpi(computer,'maci64')
% %             % Do nothing, because current_dir has been defined.
% %         end
        [layer_depth,node_depth,out_a,out_v,out_d,out_gamma,out_tau] = ...
            resp_NLH2(f_max,ppw,n_dt,boundary,N_spr,accel_in,...
            new_profile,curve,para,dir_h2rev,use_fortran);
        
        %% Process output files -- 1
        out_a = importdata('out_a.dat');
        %     out_a = fixBugOfE(out_a);
        out_v = importdata('out_v.dat');
        out_d = importdata('out_d.dat');
        out_gamma = importdata('out_gamma.dat');
        out_tau = importdata('out_tau.dat');
        layer_boundary_depth = importdata('node_depth.dat');
        layer_midpoint_depth = importdata('layer_depth.dat');
                
        layer_boundary_depth = layer_boundary_depth.';
        layer_midpoint_depth = layer_midpoint_depth.';
        
        max_a = max(abs(out_a)).';
        max_v = max(abs(out_v)).';
        max_d = max(abs(out_d)).';
        max_gamma = max(abs(out_gamma)).';
        max_tau = max(abs(out_tau)).';
        
        max_avd = [layer_boundary_depth,max_a,max_v,max_d];
        max_gt = [layer_midpoint_depth,max_gamma,max_tau];
        
        delete('out_a.dat');
        delete('out_v.dat');
        delete('out_d.dat');
        delete('layer_depth.dat');
        delete('node_depth.dat');
        delete('out_gamma.dat');
        delete('out_tau.dat');
        delete('incident.dat');
        delete('check.dat');
        delete('control.dat');
        delete('profile.dat');
        delete('curve.dat');
        delete('H2_n.dat');
        delete('t.dat');
        delete('tabk.dat');
        delete('max_d.dat');
        delete('max_v.dat');
        delete('max_tau.dat');
        delete('max_gamma.dat');
        
        if strcmpi(computer,'pcwin64')
            if use_fortran == 1
                delete('NLH2.exe');
            end
        elseif strcmpi(computer,'maci64')
            if use_fortran == 1
                delete('./NLH2');
            end
        end
        
        %         dlmwrite('max_displ_profile.dat',max_d,'delimiter','\t');
        %         dlmwrite('max_veloc_profile.dat',max_v,'delimiter','\t');
        %         dlmwrite('max_stress_profile.dat',max_tau,'delimiter','\t');
        
        cd(current_dir);  % go back to the previous "current" folder
        
        %% Process output files -- 2
        
        filename_TF_raw = sprintf('%s_nonlinear_TF_raw%s',motion_name_without_ext,ext);
        filename_TF_smoothed = sprintf('%s_nonlinear_TF_smoothed%s',motion_name_without_ext,ext);
        filename_surface_accel = sprintf('%s_accel_on_surface%s',motion_name_without_ext,ext);
        filename_new_profile = sprintf('%s_re-discretized_profile%s',motion_name_without_ext,ext);
        filename_out_a = sprintf('%s_time_history_accel%s',motion_name_without_ext,ext);
        filename_out_v = sprintf('%s_time_history_veloc%s',motion_name_without_ext,ext);
        filename_out_d = sprintf('%s_time_history_displ%s',motion_name_without_ext,ext);
        filename_out_gamma = sprintf('%s_time_history_strain%s',motion_name_without_ext,ext);
        filename_out_tau = sprintf('%s_time_history_stress%s',motion_name_without_ext,ext);
        filename_max_avd = sprintf('%s_max_a_v_d%s',motion_name_without_ext,ext);
        filename_max_gt = sprintf('%s_max_gamma_tau%s',motion_name_without_ext,ext);
        
        accel_surface = out_a(:,1);
        time_array = accel_in(:,1);
        dlmwrite(fullfile(output_dir2,filename_surface_accel),[time_array,accel_surface],'delimiter','\t','precision',7,'newline','pc');
        
        dlmwrite(fullfile(output_dir2,filename_new_profile),new_profile,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_a),out_a,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_v),out_v,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_d),out_d,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_gamma),out_gamma,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_tau),out_tau,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_max_avd),max_avd,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_max_gt),max_gt,'delimiter','\t','precision',6,'newline','pc');
        
        %% Compute non-linear transfer function
        N = length(time_array);
        dt = time_array(2) - time_array(1);
        df = 1 / (N*dt);
        freq_array = (0 : df : df*(N-1))';
        ACCEL_IN = fft(accel_in(:,2));
        ACCEL_SURFACE = fft(accel_surface(:,1));
        %     tf1 = getsmoothed(abs(ACCEL_SURFACE),df,0.3)./(2*getsmoothed(abs(ACCEL_IN),df,0.3));
        freq_array = freq_array(find(freq_array <= f_max)); % only keeps 0 -- f_max
        ACCEL_IN = ACCEL_IN(find(freq_array <= f_max));
        ACCEL_SURFACE = ACCEL_SURFACE(find(freq_array <= f_max));
        
        tf0 = abs(ACCEL_SURFACE)./(is_incident_multiplication_factor*abs(ACCEL_IN));  % raw, unsmoothed transfer function
        if use_parallel == 1
            tf2 = fasterKonnoOhmachi(abs(ACCEL_SURFACE),freq_array)...
                  ./(is_incident_multiplication_factor*fasterKonnoOhmachi(abs(ACCEL_IN),freq_array)); % smoothed transfer function
        elseif use_parallel == 0
            tf2 = fastKonnoOhmachi(abs(ACCEL_SURFACE),freq_array)...
                  ./(is_incident_multiplication_factor*fastKonnoOhmachi(abs(ACCEL_IN),freq_array)); % smoothed transfer function
        end
        
        dlmwrite(fullfile(output_dir2,filename_TF_raw),[freq_array,tf0],'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_TF_smoothed),[freq_array,tf2],'delimiter','\t','precision',6,'newline','pc');
        
        %% Fourier amplitude at ground surface
        %     FA = [freq_array(1:round(end/2)),abs(ACCEL_SURFACE(1:round(end/2)))];
        %     dlmwrite(fullfile(output_dir2,'Fourier_amplitude_at_surface.dat'),FA,'delimiter','\t');
        
        %% Export figures
        fz_axes = 12;
        fz_title = 14;
        
        % %     N = length(time_array);
        % %     dt = time_array(2) - time_array(1);
        % %     df = 1 / (N*dt);
        % %     freq_array = (0 : df : df*(N-1))';
        % %     ACCEL_IN = fft(accel_in(:,2));
        % %     ACCEL_SURFACE = fft(accel_surface(:,1));
        
        
        % transfer functions (raw and smoothed)
        fig1 = figure('visible',fig_visible);
        width = 5; height = 7;
        xLeft = (8-width)/2; yBottom = (10-height)/2;
        set(fig1,'Units','inches','paperposition',[xLeft,yBottom,width,height]);
        
        subplot(211);
        semilogx(freq_array,tf0,'color',[.6 .6 .6]); hold on;
        semilogx(freq_array,tf2,'r','linewidth',1.75);
        xlabel('Frequency (Hz)','fontsize',fz_axes);
        ylabel('Amplification Factor','fontsize',fz_axes);
        grid on;
        set(gca,'fontsize',fz_axes);
        xlim([0.1 max(freq_array)]);
        ylim([0 max(tf2(find(freq_array >= 0.1)))*1.35]);
        tf_title_text = sprintf('Nonlinear Amplification Factor\n%s',motion_name_without_ext);
        title(tf_title_text,'interpreter','none','fontsize',fz_title);
        
        subplot(212);
        plot(freq_array,tf0,'color',[.6 .6 .6]); hold on;
        plot(freq_array,tf2,'r','linewidth',1.75);
        xlabel('Frequency (Hz)','fontsize',fz_axes);
        ylabel('Amplification Factor','fontsize',fz_axes);
        grid on;
        set(gca,'fontsize',fz_axes);
        xlim([0.1 max(freq_array)]);
        ylim([0 max(tf2(find(freq_array >= 0.1)))*1.35]);
        
        tf_fig_filename = sprintf('%s_nonlinear_TF.png',motion_name_without_ext);
        saveas(fig1,fullfile(output_dir2,tf_fig_filename));
        
        % input and output accelerations
        fig2 = figure('visible',fig_visible);
        width = 5; height = 4;
        xLeft = (18.5-width)/2;  yBottom = (11-height)/2;
        set(fig2,'units','inches','paperposition',[xLeft,yBottom,width,height]);
        
        plot(time_array,accel_surface(:,1),'color',[.6 .6 .6]); hold on;
        plot(time_array,accel_in(:,2)*outcrop_conversion_factor,'r');
        ylabel('Acceleration (m/s^2)','fontsize',fz_axes);
        xlabel('Time (s)','fontsize',fz_axes);
        set(gca,'fontsize',fz_axes);
        grid on;
        xlim([0 max(time_array)]);
        accel_title_text = sprintf('Input and Output Accelerations\n%s',motion_name_without_ext);
        title(accel_title_text,'interpreter','none','fontsize',fz_title);
        
        legend('Output','Rock Outcrop','location','northeast');
        
        accel_fig_filename = sprintf('%s_accelerations.png',motion_name_without_ext);
        saveas(fig2,fullfile(output_dir2,accel_fig_filename));
        
        %% Plot max_a_v_d_gamma_tau
        fz_axes = fz_axes - 2;
        fz_title = fz_title - 1;
        
        fig3 = figure('visible',fig_visible);
        width = 12; height = 7;
        xLeft = (12.5-width)/2;  yBottom = (11-height)/2;
        set(fig3,'units','inches','paperposition',[xLeft,yBottom,width,height]);
        
        subplot(151);
        plot(max_a,layer_boundary_depth,'b.-');
        set(gca,'yDir','reverse','fontsize',fz_axes);
        ylim([0 max(layer_boundary_depth)]);
        xlabel('Max. accel. (m/s^2)','fontsize',fz_axes);
        ylabel('Depth (m)','fontsize',fz_axes);
        grid on;
        
        subplot(152);
        plot(max_v*100,layer_boundary_depth,'b.-');
        set(gca,'yDir','reverse','fontsize',fz_axes);
        ylim([0 max(layer_boundary_depth)]);
        xlabel('Max. veloc. (cm/s)');
        grid on;
        
        subplot(153);
        plot(max_d*100,layer_boundary_depth,'b.-');
        set(gca,'yDir','reverse','fontsize',fz_axes);
        ylim([0 max(layer_boundary_depth)]);
        xlabel('Max. displ. (cm)','fontsize',fz_axes);
        title(sprintf('%s,',motion_name_without_ext),'interpreter','none','fontsize',fz_title);
        grid on;
        
        subplot(154);
        plot(max_gamma*100,layer_midpoint_depth,'b.-');
        set(gca,'yDir','reverse','fontsize',fz_axes);
        ylim([0 max(layer_boundary_depth)]);
        xlabel('\gamma_{max} (%)','fontsize',fz_axes);
% %         set(get(gca,'xLabel'),'position',get(get(gca,'xLabel'),'position')+[0,4.5,0]);
        grid on;
        
        subplot(155);
        plot(max_tau/1000,layer_midpoint_depth,'b.-');
        set(gca,'yDir','reverse','fontsize',fz_axes);
        ylim([0 max(layer_boundary_depth)]);
        xlabel('\sigma_{max} (kPa)','fontsize',fz_axes);
% %         set(get(gca,'xLabel'),'position',get(get(gca,'xLabel'),'position')+[0,4.5,0]);
        grid on;
        
        max_fig_filename = sprintf('%s_max_a_v_d_gamma_tau.png',motion_name_without_ext);
        saveas(fig3,fullfile(output_dir2,max_fig_filename));
        
        % % % % %
        if strcmpi(fig_visible,'off')
            close(fig1);
            close(fig2);
            close(fig3);
        end
    end

end

end








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function calling the H2 model
%
% function [ layer_depth, node_depth, out_a, out_v, out_d, out_gamma, out_tau ] = ...
%    resp_H2(f_max, ppw, n_dt, boundary, N_spr, acc, vs_profile, curve, para, dir_H2 )
%
%       layer_depth: depth of the central point of every layer
%       node_depth:  depth of node
%       out_a:       acceleration time history of every node
%       out_v:       velocity time history of every node
%       out_d:       displacement time history of every node
%       out_gamma:   shear strain time history of every layer
%       out_tau:     shear stress time history of every layer
%       f_max:       maximum frequency in the simulation
%       ppw:         number of points per wave length
%       n_dt:        number of substep in one time step
%       boundary:    boundary condition
%       N_spr:       number of Iwan springs
%       acc:         input acceleration time history
%       vs_profile:  vs_profile data
%       curve:       modulus reduction and material damping data
%       para:        parameters for nonlinear model
%       dir_??       directory of the executable file

function [ layer_depth, node_depth, out_a, out_v, out_d, out_gamma, out_tau ] = ...
    resp_NLH2(f_max, ppw, n_dt, boundary, N_spr, acc, vs_profile, curve, para, dir_H2, use_fortran)

vs_profile(:,4) = vs_profile(:,4)/1000; % The input "profile" has a unit of rho of kg/m^3.
%                                   It is necessary to convert it to g/cm^3
%                                   to give to NLH2.

cd(dir_H2);
nlayer=size(vs_profile,1)-1;
t=acc(:,1); nt_out=length(t);
N_obs=size(curve,1); n_ma=size(curve,2)/4;

switch boundary
    case 'elastic'
        n_bound=1;
    case 'rigid'
        n_bound=2;
    case 'borehole'
        n_bound=3;
end

fid=fopen('control.dat', 'w');
fprintf(fid, '%6.1f %6.0f %6.0f %6.0f %6.0f %10.0f %6.0f %6.0f %6.0f', ...
              f_max, ppw, n_dt, n_bound, nlayer, nt_out, n_ma, N_spr, N_obs);
fclose(fid);

dlmwrite('profile.dat', vs_profile);
dlmwrite('incident.dat',acc);
dlmwrite('curve.dat', curve);
dlmwrite('H2_n.dat', para);

if strcmpi(computer,'pcwin64')
    if use_fortran == 1
        system('NLH2');
    else
        nlH2(dir_H2); % run function "nlH2" -- no input args, no output args
    end
elseif strcmpi(computer,'maci64')
    if use_fortran == 1
        system('./NLH2');
    else
        nlH2(dir_H2);
    end
end

layer_depth=load('layer_depth.dat');
node_depth=load('node_depth.dat');
out_a=load('out_a.dat');
out_v=load('out_v.dat');
out_d=load('out_d.dat');
out_gamma=load('out_gamma.dat');
out_tau=load('out_tau.dat');


end