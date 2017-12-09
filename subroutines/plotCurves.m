function plotCurves(curve_matrix)
%
% Plots G/Gmax and dampint curves
%

nr_strn_pts = size(curve_matrix,1);  % number of strain points
nr_mtrls = size(curve_matrix,2)/4; % number of different materials

strn_arrays_for_G = zeros(nr_strn_pts,nr_mtrls);
strn_arrays_for_xi = zeros(nr_strn_pts,nr_mtrls);
G_arrays = zeros(nr_strn_pts,nr_mtrls);
xi_arrays = zeros(nr_strn_pts,nr_mtrls);

for i = 1 : 1 : nr_mtrls
    strn_arrays_for_G(:,i) = curve_matrix(:,(i-1)*4+1);
    G_arrays(:,i) = curve_matrix(:,(i-1)*4+2);
    strn_arrays_for_xi(:,i) = curve_matrix(:,(i-1)*4+3);
    xi_arrays(:,i) = curve_matrix(:,(i-1)*4+4);
end

hfig = figure('unit','pixels','outerposition',[50,50,750,380]);
set(hfig,'unit','inches','paperposition',[2.25,3.5,7.5,3]);

subplot(121);
semilogx(strn_arrays_for_G,G_arrays);
xlabel('Shear Strain (%)');
ylabel('G/G_{max}');
xlim([0.0001 3]);
set(gca,'xTick',[0.0001 0.001 0.01 0.1 1 3]);
if min(min(G_arrays)) < 1.0
    ylim([min(min(G_arrays)) 1.0]);
end
grid on;

subplot(122);
semilogx(strn_arrays_for_xi,xi_arrays);
xlabel('Shear Strain (%)');
ylabel('Damping Ratio (%)');
xlim([0.0001 3]);
set(gca,'xTick',[0.0001 0.001 0.01 0.1 1 3]);
if min(min(xi_arrays)) < max(max(xi_arrays))
    ylim([min(min(xi_arrays)) max(max(xi_arrays))]);
end
grid on;
