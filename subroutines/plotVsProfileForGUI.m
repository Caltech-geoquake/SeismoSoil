function plotVsProfileForGUI(vs_profile,title_text)
% USAGE: plotVsProfileForGUI(profile_matrix,title_text)
%
% (c) Jian Shi, 7/1/2013

if nargin < 2
    title_text = ''; % empty string
end

data = vs_profile;

thickness = data(:,1);
Vs = data(:,2);
zmax = sum(thickness)+thickness(1);

plotVsProfileFromArraysWithTitle(thickness,Vs,zmax,title_text);

end

function plotVsProfileFromArraysWithTitle(h,Vs,zmax,title_text)
% USAGE: plotVsProfileFromArrays(h,Vs,zmax)
%
% (c) Jian Shi, 6/21/2012
% On 5/30/2013, function name was changed from "plotVsProfile" to
% "plotVsProfileFromArrays". (from Xiamen, China)

[x,y] = genProfilePlotArray(h,Vs,zmax);

hfig = figure('unit','pixels','outerposition',[150,150,360,600]);
set(hfig,'unit','inches','paperposition',[4.2,1.5,3.6,5]);
% % fig1 = figure; % Create a figure object named "fig1" ("handle" of the fig)
% % xSize = 3.6; ySize = 5;     %  Width = 3 in., height = 2.25 in.
% % set(fig1,'Units','inches');  %  Set units as inch
% % set(fig1,'PaperUnits','inches');  %  Set paper unit as inch
% % xLeft = (12-xSize)/2;  % xLeft: the larger the more to the right;
% % yTop = (8-ySize)/2;   % yTop: the larger the higher
% % set(fig1,'position',[xLeft yTop xSize ySize]);  
% %              % Defines where figure is shown on screen
% % set(fig1,'PaperPosition',[xLeft yTop xSize ySize]); 
% %              % Defines where figure is shown on "paper"

% h1 = axes;
plot(x,y,'k','linewidth',1.5);
set(gca, 'Ydir', 'reverse');
% set(h1, 'YAxisLocation', 'Left');
ylim([0 zmax]);
xlim([0 max(Vs)*1.2]);
xlabel('Shear-wave velocity [m/s]','fontsize',12);
ylabel('Depth [m]','fontsize',12);
set(gca,'fontsize',12);
title(title_text,'interpreter','none','fontsize',14);
grid on;

end

function [x,y] = genProfilePlotArray(h,Vs,zmax)
N = length(Vs);

z = convertThicknessToDepth(h);

x = zeros(2*N,1); % preallocation
y = zeros(2*N,1); % preallocation

for i = 1 : 1 : 2*N
    x(i) = Vs(ceil(i/2));
    if i < 2*N, y(i) = z(floor(i/2)+1); else y(i) = zmax; end
end

end

function z = convertThicknessToDepth(h)

z = zeros(length(h),1);

for i = 2 : 1 : length(h)
    z(i) = z(i-1) + h(i-1);
end

end
