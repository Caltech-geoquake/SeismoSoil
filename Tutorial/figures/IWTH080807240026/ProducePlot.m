clear all; clc; close all;

sh1 = load('IWTH080807240026.SH1.dat');
sh2 = load('IWTH080807240026.SH2.dat');

[f1,fas1] = fourierTransform(sh1,'single','abs','n');
[f2,fas2] = fourierTransform(sh2,'single','abs','n');
af = fas2./fas1;
af_smooth = fastKonnoOhmachi(af,f1);
figure;
semilogx(f1,af,'color',[.8 .8 .8]); hold on;
semilogx(f1,af_smooth,'k','linewidth',1.75);

tf_ln = load('./LN/IWTH080807240026_linear_TF.txt');
tf_eq = load('./EQ/IWTH080807240026_equivalent_linear_TF.txt');
tf_fd = load('./FD/IWTH080807240026_equivalent_linear_(FD)_TF.txt');

semilogx(tf_ln(:,1),tf_ln(:,2),'b','linewidth',1.5);
semilogx(tf_eq(:,1),tf_eq(:,2),'r','linewidth',1.5);
% semilogx(tf_fd(:,1),tf_fd(:,2),'k--','linewidth',1.5);

xlim([1e-1,50]);
ylim([0,1.1*max(tf_ln(:,2))]);

hl = legend('True (unsmoothed)','True (smoothed)','Linear',...
    'Equivalent linear');
set(hl,'fontsize',14,'location','northwest');

xlabel('Frequency (Hz)','fontsize',14);
set(gca,'fontsize',14);