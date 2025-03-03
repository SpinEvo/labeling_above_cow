clc;clear;close all;    
modmat_path = './Simulation/FigS4/';
load([modmat_path,filesep,'modmat_phase.mat']);
% Optimal CoW, RFsep = 1560us
modmat_lookup_RFsep1560 = load([modmat_path,filesep,'modmat_lookup_OptimalCoW_bipolar.mat']).modmat_lookup;
% Optimal CoW, RFsep = 1380us
modmat_lookup_RFsep1380 = load([modmat_path,filesep,'modmat_lookup_OptimalCoW_Sep1380_bipolar.mat']).modmat_lookup;
% Optimal CoW, RFsep = 1200us
modmat_lookup_RFsep1200 = load([modmat_path,filesep,'modmat_lookup_OptimalCoW_Sep1200_bipolar.mat']).modmat_lookup;

fig1 = figure;
plot(modmat_phase,modmat_lookup_RFsep1560,'LineWidth', 2);hold on;
plot(modmat_phase,modmat_lookup_RFsep1380,'LineWidth', 2);hold on;
plot(modmat_phase,modmat_lookup_RFsep1200,'LineWidth', 2);

xlabel 'Phase alternation (radians)'
ylabel 'Mz in simulation'
xticks(-pi:pi:pi); xticklabels({'-\pi','0','\pi'});xlim([-pi-0.01 pi+0.01]);
legend( 'RF = 1560μs','RF = 1380μs','RF = 1200μs' ,'location','best' ); % Legend is velocities in cm/s
axis on; box off;
set(gca,'linewidth',1.5,'FontSize',12,'FontWeight','bold');
ylim([-1 1]);
set(fig1,'Position',[0,0,500,400]); 