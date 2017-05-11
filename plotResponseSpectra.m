function plotResponseSpectra(Tn,SA,PSA,SV,PSV,SD,motion,...
                             x_axis_scale,accel_unit,motion_name)

fz_12 = 12;
fz_14 = 14;

switch accel_unit
    case 'm/s^2'
        sa_unit = accel_unit;
        sv_unit = 'm/s';
        sd_unit = 'm';
    case 'gal'
        sa_unit = accel_unit;
        sv_unit = 'cm/s';
        sd_unit = 'cm';
    case 'g'
        sa_unit = accel_unit;
        sv_unit = ' x 9.81 m/s';
        sd_unit = ' x 9.81 m';
end


fig1 = figure('unit','pixels','outerposition',[20,25,750,700]);
xSize = 9; ySize = 7;     %  Width = 3 in., height = 2.25 in.
% set(fig1,'Units','inches');  %  Set units as inch
% set(fig1,'PaperUnits','inches');  %  Set paper unit as inch
xLeft = (12-xSize)/2;  % xLeft: the larger the more to the right;
yTop = (12-ySize)/2;   % yTop: the larger the higher
% set(fig1,'position',[xLeft yTop xSize ySize]);  
%              % Defines where figure is shown on screen
set(fig1,'unit','inches','PaperPosition',[xLeft yTop xSize ySize]); 
             % Defines where figure is shown on "paper"

subplot(221);
plot(motion(:,1),motion(:,2));
xlabel('Time [sec]','fontsize',fz_12);
ylabel(sprintf('Accel. [%s]',accel_unit),'fontsize',fz_12);
xlim([min(motion(:,1)) max(motion(:,1))]);
set(gca,'fontsize',fz_12);
title(motion_name,'interpreter','none','fontsize',fz_14);
grid on;


subplot(222);
plot(Tn,SA,'linewidth',1.5);
xlabel('Period [sec]','fontsize',fz_12);
ylabel(sprintf('Spectral Acceleration [%s]',sa_unit),'fontsize',fz_12);
set(gca,'fontsize',fz_12,'xScale',x_axis_scale);
xlim([min(Tn) max(Tn)]);
if strcmpi(x_axis_scale,'log')
    new_xticklabel = get(gca,'xtick');
    set(gca,'xticklabel',new_xticklabel); % make x label non-exponential
end
grid on;

subplot(223);
plot(Tn,SV,'linewidth',1.5);
xlabel('Period [sec]','fontsize',fz_12);
ylabel(sprintf('Spectral Velocity [%s]',sv_unit),'fontsize',fz_12);
set(gca,'fontsize',fz_12,'xScale',x_axis_scale);
xlim([min(Tn) max(Tn)]);
if strcmpi(x_axis_scale,'log')
    new_xticklabel = get(gca,'xtick');
    set(gca,'xticklabel',new_xticklabel); % make x label non-exponential
end
grid on;


subplot(224);
plot(Tn,SD,'linewidth',1.5);
xlabel('Period [sec]','fontsize',fz_12);
ylabel(sprintf('Spectral Displacement [%s]',sd_unit),'fontsize',fz_12);
set(gca,'fontsize',fz_12,'xScale',x_axis_scale);
xlim([min(Tn) max(Tn)]);
if strcmpi(x_axis_scale,'log')
    new_xticklabel = get(gca,'xtick');
    set(gca,'xticklabel',new_xticklabel); % make x label non-exponential
end
grid on;