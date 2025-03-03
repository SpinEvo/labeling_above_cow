clc;clear;close all;
workingDir = './Simulation/Fig4b'; cd(workingDir); delete *asv;
utils_path = './utils_IOES'; addpath(genpath(utils_path));
modmat_path = './Simulation/Fig4b';
resultsFolder = [workingDir,filesep,'results_7Subs'];
string = {'OES_Default_Setting','IOES_Default_Setting','IOES_Optimal_Setting'};
isUnipolar = 0;
sigma = 0:0.5:6; % mm
seed_index = 10;

for seed_count = 1:100
    rng('default');
    seed_index; rng(seed_index);
    for subid = 1:7
        load(['VesLocsInfo_Sub00',num2str(subid),'.mat']);
        AngTowardsCor = VesLocsInfo.AngTowardsCor;
        AngTowardsSag = VesLocsInfo.AngTowardsSag;
        VesLocs3D = VesLocsInfo.VesLocs3D; off_res = 0;
        nVes = size(VesLocs3D,1);
        SubPath = [resultsFolder,filesep,'Sub00',num2str(subid)];
        for i_sigma = 1:length(sigma)
            distance = normrnd(0, sigma(i_sigma));
            theta = rand * 2 * pi;
            direction = [cos(theta), sin(theta)];
            perturbation = [distance * direction 0];
            perturbation = repmat(perturbation,[nVes 1]);
            VesLocs3D_shift(:,:,:,i_sigma) = VesLocs3D + perturbation;
        end
        for id_result = 1:length(string)
            VEPCASL_file_path = [SubPath,filesep,string{id_result},filesep,'VEPCASL_OES_Setup.txt'];
            CycDets = load(VEPCASL_file_path);
            % syntax: [TagAng  TagMode   vA     vB; ...]
            % Swap column order to match expected form
            CycDets = CycDets(:,[1 4 2 3]);
            for m = 1:size(VesLocs3D_shift,4)
                VesLocs3D_motion = squeeze(VesLocs3D_shift(:,:,:,m));
                [LRUnitVec, APUnitVec] = CalcObliqueUnitVecs(AngTowardsCor,AngTowardsSag);
                VesLocs = transpose(Convert3DTo2DVesLocs(VesLocs3D_motion', LRUnitVec, APUnitVec));
                if strcmp(string{id_result},'IOES_Optimal_Setting')
                    InvEffBipolarPath_0Hz = [modmat_path,filesep,'hw_PcaslParas_vs_pos_and_v_gaussian_dsv_OptimalCoW_Sep1560_bipolar_0Hz.mat'];
                else
                    InvEffBipolarPath_0Hz = [modmat_path,filesep,'hw_PcaslParas_vs_pos_and_v_gaussian_dsv_Siemens_default_3T_bipolar_0Hz.mat'];
                end
                PerfectInv = false;
                A = hw_Enc_Mtx_ProjCoW(VesLocs,CycDets,30,PerfectInv,true,isUnipolar,InvEffBipolarPath_0Hz);
                [~,E,~,~] = AssessEncoding(A); E_mean = mean(E(1:nVes));
                E_array(subid,id_result,m,seed_count) = E_mean;
            end
        end
    end
    seed_index = seed_index + 10;
end
E_array_AllSub = mean(E_array,4);
clc; rmpath(genpath(utils_path));
%% Plot
E_array_AveSub = mean(E_array_AllSub,1);
std_array = std(E_array_AllSub,1);

trans = 0.5;
fig = figure; nsub = 1;
mean_E_OES_Default = squeeze(E_array_AveSub(nsub,1,:));
mean_E_IOES_Default = squeeze(E_array_AveSub(nsub,2,:));
mean_E_IOES_Optimal = squeeze(E_array_AveSub(nsub,3,:));

std_E_OES_Default = squeeze(std_array(nsub,1,:));
std_E_IOES_Default = squeeze(std_array(nsub,2,:));
std_E_IOES_Optimal = squeeze(std_array(nsub,3,:));

scatter(sigma, mean_E_OES_Default, 30, [1 0 0], 'filled','MarkerFaceAlpha', trans); hold on;
scatter(sigma, mean_E_IOES_Default, 30, [0 0 1], 'filled','MarkerFaceAlpha', trans); hold on;
scatter(sigma, mean_E_IOES_Optimal, 30, [0.6 0.9 0.2], 'filled','MarkerFaceAlpha', trans+0.2); hold on;

plot(sigma,mean_E_OES_Default,'--','LineWidth', 2,'Color',[0 0 0 0.2]);hold on;
plot(sigma,mean_E_IOES_Default,'--','LineWidth', 2,'Color',[0 0 0 0.2]);hold on;
plot(sigma,mean_E_IOES_Optimal,'--','LineWidth', 2,'Color',[0 0 0 0.2]);hold on;

trans = 0.2;
errorbar(sigma, mean_E_OES_Default, std_E_OES_Default, 'LineStyle', 'none', 'Color', [0 0 0 trans]); hold on;
errorbar(sigma, mean_E_IOES_Default, std_E_IOES_Default, 'LineStyle', 'none', 'Color', [0 0 0 trans]); hold on;
errorbar(sigma, mean_E_IOES_Optimal, std_E_IOES_Optimal, 'LineStyle', 'none', 'Color', [0 0 0 trans]); hold on;
%legend({'OES / Default', 'IOES / Default','IOES / Optimal'},'Location','southwest','NumColumns',2);

xlabel('Magnitude of head motion (mm)');
ylabel('Mean SNR efficiency');grid on;
set(gca,'FontWeight','bold','FontSize',12,'LineWidth',1.5)
ylim([0.3 0.8]);
set(fig,'Position',[0,0,500,250]);