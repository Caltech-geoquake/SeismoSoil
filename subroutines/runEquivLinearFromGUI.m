function [ok_to_proceed,h_running] = runEquivLinearFromGUI(vs_profile,curve,nr_motion,motion,motion_name,output_dir,...
        factor_rho,factor_xi,unit_factor_accel,bedrock_type,motion_type,fig_visible_option,use_parallel)
%
h_running = 0;

%% Unit conversions -- soil profile
rho = vs_profile(:,4);
xi = vs_profile(:,3);
rho = rho ./ factor_rho;
xi = xi ./ factor_xi;
vs_profile(:,4) = rho;
vs_profile(:,3) = xi;

if use_parallel == 1
    nr_cores = inf; % using a maximum of "nr_cores" workers or threads
elseif use_parallel == 0
    nr_cores = 0; % do not use multiple cores/threads
end

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
    
    parfor (ii = 1 : nr_motion, nr_cores)
        fprintf('Ground motion No. %d\n',ii);
        
        current_motion = motion{ii};
        [~,motion_name_without_ext,ext] = fileparts(motion_name{ii});
        mkdir(output_dir,motion_name_without_ext);
        output_dir2 = fullfile(output_dir,motion_name_without_ext);
        
        accel_temp = current_motion(:,2);
        accel_temp = accel_temp ./ unit_factor_accel;
        current_motion(:,2) = accel_temp;
        accel_incident = current_motion;
        accel_incident(:,2) = accel_incident(:,2)/accel_division_factor;
        
        output_or_not = 'n';
        [freq_array,tf,t_out,accel_on_surface,...
            new_profile,out_a,out_v,out_d,out_gamma,out_tau,max_avd,max_gt] ...
            = equivLinSiteRespVec(vs_profile,accel_incident,curve,'off',boundary,output_or_not);
        
        filename_TF = sprintf('%s_equivalent_linear_TF%s',motion_name_without_ext,ext);
        filename_surface_accel = sprintf('%s_accel_on_surface%s',motion_name_without_ext,ext);
        filename_new_profile = sprintf('%s_re-discretized_profile%s',motion_name_without_ext,ext);
        filename_out_a = sprintf('%s_time_history_accel%s',motion_name_without_ext,ext);
        filename_out_v = sprintf('%s_time_history_veloc%s',motion_name_without_ext,ext);
        filename_out_d = sprintf('%s_time_history_displ%s',motion_name_without_ext,ext);
        filename_out_gamma = sprintf('%s_time_history_strain%s',motion_name_without_ext,ext);
        filename_out_tau = sprintf('%s_time_history_stress%s',motion_name_without_ext,ext);
        filename_max_avd = sprintf('%s_max_a_v_d%s',motion_name_without_ext,ext);
        filename_max_gt = sprintf('%s_max_gamma_tau%s',motion_name_without_ext,ext);
        
        dlmwrite(fullfile(output_dir2,filename_TF),[freq_array,tf],'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_surface_accel),[t_out,accel_on_surface],'delimiter','\t','precision',7,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_new_profile),new_profile,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_a),out_a,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_v),out_v,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_d),out_d,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_gamma),out_gamma,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_out_tau),out_tau,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_max_avd),max_avd,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite(fullfile(output_dir2,filename_max_gt),max_gt,'delimiter','\t','precision',6,'newline','pc');
        
        %% Plot transfer function
        fz_axes = 12;
        fz_title = 14;
        
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
        tf_title_text = sprintf('Equivalent Linear Amplification Factor\n%s',motion_name_without_ext);
        title(tf_title_text,'interpreter','none','fontsize',fz_title);
        
        subplot(212);
        plot(freq_array,tf,'k','linewidth',1.5);
        xlabel('Frequency (Hz)','fontsize',fz_axes);
        ylabel('Amplification Factor','fontsize',fz_axes);
        grid on;
        set(gca,'fontsize',fz_axes);
        xlim([min(freq_array) max(freq_array)]);
        
        tf_fig_filename = sprintf('%s_equivalent_linear_TF.png',motion_name_without_ext);
        saveas(fig1,fullfile(output_dir2,tf_fig_filename));
        
        %% Plot accelerations
        fig2 = figure('visible',fig_visible_option); % Create a figure object named "fig1" ("handle" of the fig)
        width = 5; height = 4;
        xLeft = (18.5-width)/2;  yBottom = (11-height)/2;
        set(fig2,'units','inches','paperposition',[xLeft,yBottom,width,height]);
        
        if max(abs(accel_on_surface)) > max(abs(accel_incident(:,2)*outcrop_conversion_factor))
            plot(t_out,accel_on_surface,'color',[255,127,14]/255,'linewidth',1.75); hold on;
            plot(accel_incident(:,1),accel_incident(:,2)*outcrop_conversion_factor,'color',[31,119,180]/255,'linewidth',1.75);
            h_legend = legend('Output','Rock outcrop','location','northeast');%,'fontsize',fz_axes);
        else
            plot(accel_incident(:,1),accel_incident(:,2)*outcrop_conversion_factor,'color',[31,119,180]/255,'linewidth',1.75); hold on;
            plot(t_out,accel_on_surface,'color',[255,127,14]/255,'linewidth',1.75);
            h_legend = legend('Rock outcrop','Output','location','northeast');%,'fontsize',fz_axes);
        end
        ylabel('Acceleration (m/s^2)','fontsize',fz_axes);
        xlabel('Time (s)','fontsize',fz_axes);
        set(gca,'fontsize',fz_axes);
        xlim([0 max(t_out)]);
        grid on;
        
        set(h_legend,'fontsize',fz_axes); % actually the fontsize of legend is by default the same as axes
        accel_title_text = sprintf('Input and Output Accelerations\n%s',motion_name_without_ext);
        title(accel_title_text,'interpreter','none','fontsize',fz_title);
        
        accel_fig_filename = sprintf('%s_accelerations.png',motion_name_without_ext);
        saveas(fig2,fullfile(output_dir2,accel_fig_filename));
        
        %% Plot max_a_v_d_gamma_tau
        fz_axes = fz_axes - 2;
        fz_title = fz_title - 1;
        
        layer_boundary_depth = max_avd(:,1); % nr_layer x 1
        max_a = max_avd(:,2);
        max_v = max_avd(:,3);
        max_d = max_avd(:,4);
        layer_midpoint_depth = max_gt(:,1);  % (nr_layer-1) x 1
        max_gamma = max_gt(:,2);
        max_tau = max_gt(:,3);
        
        fig3 = figure('visible',fig_visible_option);
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
        title(motion_name_without_ext,'interpreter','none','fontsize',fz_title);
        grid on;
        
        subplot(154);
        plot(max_gamma*100,layer_midpoint_depth,'b.-');
        set(gca,'yDir','reverse','fontsize',fz_axes);
        ylim([0 max(layer_boundary_depth)]);
        xlabel('\gamma_{max} (%)','fontsize',fz_axes);
        % set(get(gca,'xLabel'),'position',get(get(gca,'xLabel'),'position')+[0,4.5,0]);
        grid on;
        
        subplot(155);
        plot(max_tau/1000,layer_midpoint_depth,'b.-');
        set(gca,'yDir','reverse','fontsize',fz_axes);
        ylim([0 max(layer_boundary_depth)]);
        xlabel('\sigma_{max} (kPa)','fontsize',fz_axes);
        % set(get(gca,'xLabel'),'position',get(get(gca,'xLabel'),'position')+[0,4.5,0]);
        grid on;
        
        max_fig_filename = sprintf('%s_max_a_v_d_gamma_tau.png',motion_name_without_ext);
        saveas(fig3,fullfile(output_dir2,max_fig_filename));
        
        %%
        if strcmpi(fig_visible_option,'off')
            close(fig1);
            close(fig2);
            close(fig3);
        end
        
        fz_axes = fz_axes + 2;
        fz_title = fz_title + 1;
    end
end
