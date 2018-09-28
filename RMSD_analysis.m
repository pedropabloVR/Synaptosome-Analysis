% This script reads in cropped localisation files (from synaptoCrop, in
% rapidstorm format) and estimate the radius of clusters using simple RMSD
% calculation for all individual regions. The RMSDs for all regions are
% plotted and a mean ± SD is returned.
% 
% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated on 19 Sept 2018


%% Parameters

clear all
close all
clc

format = 'thunderstorm'; % reconstruction software used

% Path to folders of the three conditions (which all contain three
% subfolders: 'RC','GC' and 'BC' as written away by synaptoCrop)

% dir_phys  = 'F:\Data\Synaptosomes\Experiment_4C\Results\Results_combined\ripley\phys';
% dir_egta  = 'F:\Data\Synaptosomes\Experiment_4C\Results\Results_combined\ripley\egta';
% dir_egtak = 'F:\Data\Synaptosomes\Experiment_4C\Results\Results_combined\ripley\egtak';

dir_phys  = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_combined/ripley/phys';
dir_egta  = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_combined/ripley/egta';
dir_egtak = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_combined/ripley/egtak';

% dir_phys  = fullfile(pwd,'testdata/condition1_phys');  % path to testdata condition 1
% dir_egta  = fullfile(pwd,'testdata/condition2_egta');  % path to testdata condition 2
% dir_egtak = fullfile(pwd,'testdata/condition3_egtak'); % path to testdata condition 3

% Results will be saved in a new folder 'rmsd' in this directory
%output_dir = 'F:\Data\Synaptosomes\Experiment_4C\Results\Results_combined';
output_dir = '/Volumes/WD Ezra/Dump';

flagsave = 1; % 1 to write away results

stats = 0; % 1 to also do statistical analysis (One-way ANOVA)
% Note that the conclusions from the statistical analysis are only
% meaningful if the conditions for ANOVA are met (data within each group
% follow a normal distribution and their variances are the same).

%% Creating output folder in output_dir

% Create new output folder
path_output = fullfile(output_dir,'rmsd');
if exist(path_output, 'dir')
    opts.Interpreter = 'tex';
    opts.Default = 'Continue';
    quest = '\fontsize{12}An output folder ''rmsd'' already exists. If you continue, data in this folder might be overwritten.';
    answer = questdlg(quest,'Message','Cancel','Continue',opts);
    if strcmp(answer,'Continue')
        mkdir(path_output);
    else
        return
    end
else
    mkdir(path_output);
end


%% Get RMSD

tic

% Get rmsd for PHYS
[rmsd_phys_RC,~,~]  = calculateRMSDs(fullfile(dir_phys, 'RC'),format);
[rmsd_phys_GC,~,~]  = calculateRMSDs(fullfile(dir_phys, 'GC'),format);
[rmsd_phys_BC,~,~]  = calculateRMSDs(fullfile(dir_phys, 'BC'),format);

% Get rmsd for EGTA
[rmsd_egta_RC,~,~]  = calculateRMSDs(fullfile(dir_egta, 'RC'),format);
[rmsd_egta_GC,~,~]  = calculateRMSDs(fullfile(dir_egta, 'GC'),format);
[rmsd_egta_BC,~,~]  = calculateRMSDs(fullfile(dir_egta, 'BC'),format);

% Get rmsd for EGTA/K+
[rmsd_egtak_RC,~,~] = calculateRMSDs(fullfile(dir_egtak,'RC'),format);
[rmsd_egtak_GC,~,~] = calculateRMSDs(fullfile(dir_egtak,'GC'),format);
[rmsd_egtak_BC,~,~] = calculateRMSDs(fullfile(dir_egtak,'BC'),format);

if flagsave
    path_rmsd = fullfile(path_output,'rmsd.mat');
    save(path_rmsd,'rmsd_phys_RC','rmsd_phys_GC','rmsd_phys_BC',...
                   'rmsd_egta_RC','rmsd_egta_GC','rmsd_egta_BC',...
                   'rmsd_egtak_RC','rmsd_egtak_GC','rmsd_egtak_BC');
end


%% Save restults as tab-delimited txt file (to read in in R) 

clear var phys egta egtak
[phys{1:size(rmsd_phys_RC,2)}] = deal('phys');
[egta{1:size(rmsd_egta_RC,2)}] = deal('egta');
[egtak{1:size(rmsd_egtak_RC,2)}] = deal('egtak');
column_conditions = [phys egta egtak];

results_rmsd = [[rmsd_phys_RC' rmsd_phys_GC' rmsd_phys_BC'];
                [rmsd_egta_RC' rmsd_egta_GC' rmsd_egta_BC'];
                [rmsd_egtak_RC' rmsd_egtak_GC' rmsd_egtak_BC']];
results_rmsd = array2table(results_rmsd,'VariableNames',{'rmsd_RC','rmsd_GC','rmsd_BC'});
results_rmsd.condition = column_conditions';
results_rmsd = results_rmsd(:,[end 1:end-1]);

path_results_rmsd = fullfile(path_output,'results_rmsd.txt');
writetable(results_rmsd,path_results_rmsd,'Delimiter','\t');


%% Plot results red channel (mCling)

clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(rmsd_phys_RC,2)}]   = deal('PHYS');
[condition_EGTA{1:size(rmsd_egta_RC,2)}]   = deal('EGTA');
[condition_EGTAK{1:size(rmsd_egtak_RC,2)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

fig1 = figure;
subplot(121) % beeswarm plot
h = plotSpread({rmsd_phys_RC,rmsd_egta_RC,rmsd_egtak_RC},'xNames',{'PHYS','EGTA','EGTA/K+'});
set(h{1},'color','k','markersize',10);
ylabel('RMSD (nm)');
title('RMSD mCling');
set(gca,'fontsize',14);

subplot(122) % boxplot
boxplot([rmsd_phys_RC rmsd_egta_RC rmsd_egtak_RC],conditions);
ylabel('RMSD (nm)');
title('RMSD mCling');
set(gca,'fontsize',14);

if flagsave; savefig(fig1,fullfile(path_output,'rmsd_mcling.fig')); end


%% Plot results green channel (a-synuclein)

clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(rmsd_phys_GC,2)}]   = deal('PHYS');
[condition_EGTA{1:size(rmsd_egta_GC,2)}]   = deal('EGTA');
[condition_EGTAK{1:size(rmsd_egtak_GC,2)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

fig2 = figure;
subplot(121) % beeswarm plot
h = plotSpread({rmsd_phys_GC,rmsd_egta_GC,rmsd_egtak_GC},'xNames',{'PHYS','EGTA','EGTA/K+'});
set(h{1},'color','k','markersize',10);
ylabel('RMSD (nm)');
title('RMSD a-synuclein');
set(gca,'fontsize',14);

subplot(122) % boxplot
boxplot([rmsd_phys_GC rmsd_egta_GC rmsd_egtak_GC],conditions);
ylabel('RMSD (nm)');
title('RMSD a-synuclein');
set(gca,'fontsize',14);

if flagsave; savefig(fig2,fullfile(path_output,'rmsd_a-synuclein.fig')); end


%% Plot results blue channel (VAMP2)

clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(rmsd_phys_BC,2)}]   = deal('PHYS');
[condition_EGTA{1:size(rmsd_egta_BC,2)}]   = deal('EGTA');
[condition_EGTAK{1:size(rmsd_egtak_BC,2)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

fig3 = figure;
subplot(121) % beeswarm plot
h = plotSpread({rmsd_phys_BC,rmsd_egta_BC,rmsd_egtak_BC},'xNames',{'PHYS','EGTA','EGTA/K+'});
set(h{1},'color','k','markersize',10);
ylabel('RMSD (nm)');
title('RMSD VAMP2');
set(gca,'fontsize',14);

subplot(122) % boxplot
boxplot([rmsd_phys_BC rmsd_egta_BC rmsd_egtak_BC],conditions);
ylabel('RMSD (nm)');
title('RMSD VAMP2');
set(gca,'fontsize',14);

if flagsave; savefig(fig3,fullfile(path_output,'rmsd_VAMP2.fig')); end

toc