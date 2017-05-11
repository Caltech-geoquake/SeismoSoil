function plotFourierSpectrumAndOrSave(motion,x_axis_scale,y_axis_scale,side_option,...
            complex_option,smooth_option,motion_filename,save_option,...
            grid_on_option,motion_dir)

[f_array,spectrum] = fourierTransform(motion,side_option,complex_option);
if strcmpi(smooth_option,'regular')
    df = f_array(2) - f_array(1);
    smoothed_spectrum = getsmoothed(abs(spectrum),df,0.3);
else if strcmpi(smooth_option,'konnoOhmachi')
    smoothed_spectrum = fasterKonnoOhmachi(abs(spectrum),f_array);
    end
end
        
figure;
if strcmpi(smooth_option,'none') ~= 1  % if regular smooth or konno-ohmachi smooth is chosen
    plot(f_array,abs(spectrum),'color',[.3 .3 .3]);
    set(gca,'xScale',x_axis_scale);
    set(gca,'yScale',y_axis_scale);
    hold on;
    plot(f_array,smoothed_spectrum,'r','linewidth',1.75);
    legend('Original','Smoothed','location','northwest');
else
    plot(f_array,abs(spectrum),'k');
    set(gca,'xScale',x_axis_scale);
    set(gca,'yScale',y_axis_scale);
end
xlim([min(f_array) max(f_array)]);
set(gca,'fontsize',12);
if strcmpi(x_axis_scale,'log')
    new_xticklabel = get(gca,'xtick');
    set(gca,'xticklabel',new_xticklabel); % make x label non-exponential
end
xlabel('Frequency (Hz)','fontsize',12);
ylabel('Amplitude','fontsize',12);
if grid_on_option == 1
    grid on;
end
title(motion_filename,'interpreter','none','fontsize',14);

if save_option == 1
    [f_array,spectrum_cmp] = fourierTransform(motion,side_option,'complex');
    if strcmpi(smooth_option,'none') ~= 1  % if regular smooth or konno-ohmachi smooth is chosen
        spectrum_abs = smoothed_spectrum;
    else
        [f_array,spectrum_abs] = fourierTransform(motion,side_option,'abs');
    end
        
    [a0_tmp,motion_filename_without_ext,ext_tmp] = fileparts(motion_filename);
    fname_str1 = sprintf('%s_Fourier_spectrum_absolute%s',motion_filename_without_ext,ext_tmp);
    fname_str2 = sprintf('%s_Fourier_spectrum_complex%s',motion_filename_without_ext,ext_tmp);
    
    dlmwrite(fullfile(motion_dir,fname_str1),[f_array,spectrum_abs],'delimiter','\t','precision',6);
    dlmwrite(fullfile(motion_dir,fname_str2),[f_array,spectrum_cmp],'delimiter','\t','precision',6);
end

end