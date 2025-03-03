clc;clear;close all;
workingDir = 'H:\Research-2024\ZIC\Zhensen\ZIC\ProjCoW\code_sharing\prep_scan'; cd(workingDir);
utils_path = 'H:\Research-2024\ZIC\Zhensen\ZIC\ProjCoW\code_sharing\utils_IOES'; addpath(genpath(utils_path));
ASL_Setting = {'Default_Setting','Optimal_Setting'};

Grad_Mode = {'Bipolar','Unipolar'};
Method = {'OES','IOES'};

B0 = 3; % Scanner: 3 Tesla / 7 Tesla

%% CoW coordinates
%##########################################################################
AngTowardsCor = -22;
AngTowardsSag = 0;
% H-positive, F-negative
% L-positive, R-negative
% A-positive, P-negative
VesLocs3D = [-8 53 24; % RACA
            -51 31 15; % RMCA-1
            -52 19 10; % RMCA-2
            -51 11 7;  % RMCA-3
            40 33 16; % LMCA-1
            46 26 13; % LMCA-2
            49 17 10; % LMCA-3
           -13 -13 -3;  % RPCA
            12 -10 -1; % LPCA
            -5 62 28]; % LACA

%off_res = [200 300 212 230 189 100 202 100 280]/4096*pi; % from GRE field map
off_res = 0;
%##########################################################################
vesloc = VesLocs3D';
[LRUnitVec, APUnitVec] = CalcObliqueUnitVecs(AngTowardsCor,AngTowardsSag);
vesloc = Convert3DTo2DVesLocs(vesloc, LRUnitVec, APUnitVec);

%% OES - Tom's version
for id_ASL = 1:length(ASL_Setting)
    if strcmp(ASL_Setting{id_ASL},'Default_Setting')
        modmat_path = [utils_path,filesep,ASL_Setting{id_ASL},filesep,Grad_Mode{1}];
        load([modmat_path,filesep,'modmat_lookup.mat']);
        load([modmat_path,filesep,'modmat_phase.mat']);
    elseif strcmp(ASL_Setting{id_ASL},'Optimal_Setting')
        modmat_path = [utils_path,filesep,ASL_Setting{id_ASL},filesep,Grad_Mode{1}];
        load([modmat_path,filesep,'modmat_lookup.mat']);
        load([modmat_path,filesep,'modmat_phase.mat']);
    end
    
    for id_method = 1:length(Method)
        if strcmp(Method{id_method},'OES') && strcmp(ASL_Setting{id_ASL},'Default_Setting')
            tic;
            [cost,condition,min_half_lambda,ang_locs_tagmode] = fun_CalEncodCost(vesloc,modmat_phase,modmat_lookup,off_res);
            elapsedTime = toc;
            fprintf(['\n',Method{id_method},'_',ASL_Setting{id_ASL},' Elapsed time: %.2f seconds\n'], elapsedTime);

        elseif strcmp(Method{id_method},'IOES')
            niters = 100; tic;
            cost = fun_CalEncodCost(vesloc,modmat_phase,modmat_lookup,off_res);
            [~,condition,min_half_lambda,ang_locs_tagmode] = fun_FindOptimalEncod(vesloc,modmat_phase,modmat_lookup,niters,cost,off_res);
            ang_locs_tagmode(end-1,:) = [0 -100 100 0]; % tag all
            ang_locs_tagmode(end,:) = [0 -100 100 1]; % control all
            elapsedTime = toc;
            fprintf(['\n',Method{id_method},'_',ASL_Setting{id_ASL},' Elapsed time: %.2f seconds\n'], elapsedTime);
        elseif strcmp(Method{id_method},'OES') && strcmp(ASL_Setting{id_ASL},'Optimal_Setting')
            continue;
        end
        if mean(off_res)==0 && B0==3
            ang_locs_tagmode(end-1,:) = [0 -100 100 0]; % tag all
            ang_locs_tagmode(end,:) = [0 -100 100 1]; % control all
        end
        if strcmp(Method{id_method},'OES')
            save_path = [workingDir,filesep,'results',filesep,'OES_Default_Optimal_Setting'];
        else
            save_path = [workingDir,filesep,'results',filesep,Method{id_method},'_',ASL_Setting{id_ASL}];
        end
        if exist(save_path,'dir')==7
            rmdir(save_path,'s');
        end
        mkdir(save_path);

        writematrix(ang_locs_tagmode,[save_path,filesep,'VEPCASL_OES_Setup.txt']);

    end
end
rmpath(genpath(utils_path));