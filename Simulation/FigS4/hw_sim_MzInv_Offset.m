clc;clear;close all;
tom_git_path = './bloch_sim-main';
addpath(genpath(tom_git_path));
isUnipolar = 0;
v = 30/100; % 30cm/s = 0.3m/s
ORFreq = [0 50 100]; % Hz
%% Simulation
% Optimal CoW, RFsep = 1560us
Gmax = 9; Gmean = 0.45; FA = 30;
RF_dur = 870*10^-6; RF_sep = 1560*10^-6;
[IE_1560,ph_alt] = hw_fun_sim_offset(v,ORFreq,Gmax,Gmean,FA,RF_dur,RF_sep,isUnipolar);

% Optimal CoW, RFsep = 1380us
Gmax = 9; Gmean = 0.45; FA = 30;
RF_dur = 870*10^-6; RF_sep = 1380*10^-6;
[IE_1380,~] = hw_fun_sim_offset(v,ORFreq,Gmax,Gmean,FA,RF_dur,RF_sep,isUnipolar);

% Optimal CoW, RFsep = 1200us
Gmax = 9; Gmean = 0.45; FA = 30;
RF_dur = 870*10^-6; RF_sep = 1200*10^-6;
[IE_1200,~] = hw_fun_sim_offset(v,ORFreq,Gmax,Gmean,FA,RF_dur,RF_sep,isUnipolar);

IE = cat(3,IE_1560,IE_1380,IE_1200);
%% Plot
clc;close all;
fig1 = figure;
for i = 1:size(IE,3)
    plot(ph_alt,squeeze(IE(1,:,i)),'linewidth',2.5);hold on;
end
xlabel 'Phase alternation (radians)'
ylabel 'Inversion Efficiency'
xticks(-2*pi:pi:2*pi); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});xlim([-2*pi-0.01 2*pi+0.01]);ylim([0 1]);
%legend('RF = 1560μs','RF = 1380μs','RF = 1200μs','location','best','FontSize',8,'FontWeight','bold');
axis on; box off; title('No off-resonance','FontSize',12,'FontWeight','bold');
set(gca,'linewidth',1.5,'FontSize',12,'FontWeight','bold');
set(fig1,'Position',[0,0,500,300]); 

fig2 = figure;
for i = 1:size(IE,3)
    plot(ph_alt,squeeze(IE(2,:,i)),'linewidth',2.5);hold on;
end
xlabel 'Phase alternation (radians)'
ylabel 'Inversion Efficiency'
xticks(-2*pi:pi:2*pi); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});xlim([-2*pi-0.01 2*pi+0.01]);ylim([0 1]);
%legend('RF = 1560μs','RF = 1380μs','RF = 1200μs','location','best','FontSize',8,'FontWeight','bold');
axis on; box off; title('Off-resonance: 50 Hz','FontSize',12,'FontWeight','bold');
set(gca,'linewidth',1.5,'FontSize',12,'FontWeight','bold');
set(fig2,'Position',[0,0,500,300]); 

fig3 = figure;
for i = 1:size(IE,3)
    plot(ph_alt,squeeze(IE(3,:,i)),'linewidth',2.5);hold on;
end
xlabel 'Phase alternation (radians)'
ylabel 'Inversion Efficiency'
xticks(-2*pi:pi:2*pi); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});xlim([-2*pi-0.01 2*pi+0.01]);ylim([0 1]);
%legend('RF = 1560μs','RF = 1380μs','RF = 1200μs','location','best','FontSize',8,'FontWeight','bold');
axis on; box off; title('Off-resonance: 100 Hz','FontSize',12,'FontWeight','bold');
set(gca,'linewidth',1.5,'FontSize',12,'FontWeight','bold');
set(fig3,'Position',[0,0,500,300]); 

rmpath(genpath(tom_git_path));

%% Plot together
clc;close all;
fig4 = figure;COLOR = [1 0 0;0 0 1;0.4 0.7 0.5]; transp = 0.7;
for i = 1:size(IE,3)
    plot(ph_alt,squeeze(IE(1,:,i)),'Color',[COLOR(i,:),transp],'linewidth',2.5);hold on;
    plot(ph_alt,squeeze(IE(2,:,i)),'--','Color',[COLOR(i,:),transp],'linewidth',2.5);hold on;
    plot(ph_alt,squeeze(IE(3,:,i)),':','Color',[COLOR(i,:),transp],'linewidth',2.5);hold on;
end
xlabel 'Phase alternation (radians)'
ylabel 'Inversion Efficiency'
xticks(-2*pi:pi:2*pi); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});xlim([-2*pi-0.01 2*pi+0.01]);ylim([0 1]);
%legend('RF = 1560μs','RF = 1380μs','RF = 1200μs','location','best','FontSize',8,'FontWeight','bold');
axis on; box off;
set(gca,'linewidth',1.5,'FontSize',12,'FontWeight','bold');
set(fig4,'Position',[0,0,500,200]); 

xticks(-2*pi:pi:2*pi); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});xlim([-2*pi-0.01 2*pi+0.01]);ylim([0 1]);
