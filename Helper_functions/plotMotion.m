function [hfig,hlines] = plotMotion(accel,unit,unit_conversion_factor,title_text,show_RMS_accel)
% [hfig,hlines] = plotMotion(accel,unit,unit_conversion_factor,title_text,show_RMS_accel)
%  
%  INPUTS:
%     accel: acceleration; a matrix in two columns
%     unit: a string, for example, "m/s/s", or "gal"
%     unit_conversion_factor: 1 if the unit of accel is m/s/s;
%                             100 for cm/s/s; 1/9.81 for g
%     title_text: a string that goes on top of the figure as the title
%     show_RMS_accel: 0, don't show; 1, show RMS acceleration as two
%                     horizontal lines in the first subfig
%
%  (c) Jian Shi
%  Original creation date unknown. Modified 11/26/2015 (Thanksgiving Day)

if nargin < 5
    show_RMS_accel = 0;
end
if nargin < 4
    title_text = 'Accelerogram';
end
if nargin < 3
    unit_conversion_factor = 1;
end
if nargin < 2
    unit = 'm/s/s';
end

t = accel(:,1);
a = accel(:,2);

factor = unit_conversion_factor;

%% Find PGA and the corresponding index
abs_a = abs(a);
[PGA,pga_index] = max(abs_a);

%% Calculate RMS acceleration
if show_RMS_accel == 1
    dt = t(2)-t(1);
    a_sq = a.^2;
    Td = t(end);
    a_rms = sqrt(1/Td * sum(a_sq)*dt);
end

%% Numerical integration to get velocity and displacement
% a = a - a(1);
[v,u] = numInt([t,a]);

%% Plot figures
hfig = figure('unit','pixels','outerposition',[100,25,500,715]);
set(hfig,'unit','inches','paperposition',[2.5,1,5,7]);
% % xSize = 5; ySize = 7;     %  Width = 3 in., height = 2.25 in.
% % set(fig1,'Units','inches');  %  Set units as inch
% % set(fig1,'PaperUnits','inches');  %  Set paper unit as inch
% % xLeft = (10-xSize)/2;  % xLeft: the larger the more to the right;
% % yBottom = (7-ySize)/2;   % yTop: the larger the higher
% % set(hfig,'position',[xLeft yBottom xSize ySize]); % Defines where figure is shown on screen
% % set(hfig,'PaperPosition',[xLeft yBottom xSize ySize]); % Defines where figure is shown on "paper"

subplot(311);
if show_RMS_accel == 1
    plot(t,a_rms*ones(size(t)),'r--'); hold on;
    plot(t,-a_rms*ones(size(t)),'r--'); hold on;
end
h1 = plot(t,a,'b'); hold on;
plot(t(pga_index),a(pga_index),'ro','linewidth',1.75);
grid on;
% annotation('textbox',[.3 .8 .1 .1],'String',sprintf('PGA = %.2f (unit)',PGA),'backgroundcolor','w');
text(t(min(pga_index+round(length(t)/40),length(t))),a(pga_index),sprintf('PGA = %.2f %s',PGA,unit),'fontsize',12,'backgroundcolor','w');
xlim([min(t) max(t)]);
xlabel('Time (s)');
ylabel(sprintf('Acceleration [%s]',unit));
title(title_text,'fontsize',14,'interpreter','none');

subplot(312);
h2 = plot(t,v(:,2)/factor*100,'b'); % division by "factor" makes the unit into meter
grid on;
xlim([min(t) max(t)]);
xlabel('Time [s]');
ylabel('Velocity [cm/s]');

subplot(313);
h3 = plot(t,u(:,2)/factor*100,'b');
grid on;
xlim([min(t) max(t)]);
xlabel('Time [s]');
ylabel('Displacement [cm]');

hlines = {h1;h2;h3};

end