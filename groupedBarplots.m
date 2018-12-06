
clear all
close all
clc


% Load data
% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/results_combined_after_overlap_threshold.mat');
% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_combined/ripley/clustersizes.mat');
% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_combined/ripley/interclusterdist.mat');

load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/results_combined_after_overlap_threshold.mat');
load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/Results_combined/ripley/clustersizes.mat');
load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/Results_combined/ripley/interclusterdist.mat');

output_dir = '/Users/Ezra/Desktop';

temperature = '4';

%% Creating output folder in output_dir

% Create new output folder
path_output = fullfile(output_dir,'grouped_barplots');
if exist(path_output, 'dir')
    opts.Interpreter = 'tex';
    opts.Default = 'Continue';
    quest = '\fontsize{12}An output folder ''Grouped barplots'' already exists. If you continue, data in this folder might be overwritten.';
    answer = questdlg(quest,'Message','Cancel','Continue',opts);
    if strcmp(answer,'Continue')
        mkdir(path_output);
    else
        return
    end
else
    mkdir(path_output);
end

%% Plot

T = results_combined_after_filtering;

% Get mean and sd of all the data to be included in grouped barplots
[area_mCling_mean, area_mCling_std] = ...
    getMeanAndStdConditions(T.Area,T.condition,{'phys','egta','egtak'});

[overlapWithGreen_mean, overlapWithGreen_std] = ...
    getMeanAndStdConditions(T.OverlapWithGreen,T.condition,{'phys','egta','egtak'});

[overlapWithBlue_mean, overlapWithBlue_std] = ...
    getMeanAndStdConditions(T.OverlapWithBlue,T.condition,{'phys','egta','egtak'});

[weightedOverlapWithGreen_mean, weightedOverlapWithGreen_std] = ...
    getMeanAndStdConditions(T.WeightedOverlapWithGreen,T.condition,{'phys','egta','egtak'});

[weightedOverlapWithBlue_mean, weightedOverlapWithBlue_std] = ...
    getMeanAndStdConditions(T.WeightedOverlapWithBlue,T.condition,{'phys','egta','egtak'});

size_phys_RC_mean = mean(clustersize_phys_RC);
size_phys_GC_mean = mean(clustersize_phys_GC);
size_phys_BC_mean = mean(clustersize_phys_BC);
interclusterdist_phys_RG_mean = mean(interclusterdist_phys_RG);
interclusterdist_phys_RB_mean = mean(interclusterdist_phys_RB);
interclusterdist_phys_GB_mean = mean(interclusterdist_phys_GB);
size_phys_RC_std = std(clustersize_phys_RC);
size_phys_GC_std = std(clustersize_phys_GC);
size_phys_BC_std = std(clustersize_phys_BC);
interclusterdist_phys_RG_std = std(interclusterdist_phys_RG);
interclusterdist_phys_RB_std = std(interclusterdist_phys_RB);
interclusterdist_phys_GB_std = std(interclusterdist_phys_GB);

size_egta_RC_mean = mean(clustersize_egta_RC);
size_egta_GC_mean = mean(clustersize_egta_GC);
size_egta_BC_mean = mean(clustersize_egta_BC);
interclusterdist_egta_RG_mean = mean(interclusterdist_egta_RG);
interclusterdist_egta_RB_mean = mean(interclusterdist_egta_RB);
interclusterdist_egta_GB_mean = mean(interclusterdist_egta_GB);
size_egta_RC_std = std(clustersize_egta_RC);
size_egta_GC_std = std(clustersize_egta_GC);
size_egta_BC_std = std(clustersize_egta_BC);
interclusterdist_egta_RG_std = std(interclusterdist_egta_RG);
interclusterdist_egta_RB_std = std(interclusterdist_egta_RB);
interclusterdist_egta_GB_std = std(interclusterdist_egta_GB);

size_egtak_RC_mean = mean(clustersize_egtak_RC);
size_egtak_GC_mean = mean(clustersize_egtak_GC);
size_egtak_BC_mean = mean(clustersize_egtak_BC);
interclusterdist_egtak_RG_mean = mean(interclusterdist_egtak_RG);
interclusterdist_egtak_RB_mean = mean(interclusterdist_egtak_RB);
interclusterdist_egtak_GB_mean = mean(interclusterdist_egtak_GB);
size_egtak_RC_std = std(clustersize_egtak_RC);
size_egtak_GC_std = std(clustersize_egtak_GC);
size_egtak_BC_std = std(clustersize_egtak_BC);
interclusterdist_egtak_RG_std = std(interclusterdist_egtak_RG);
interclusterdist_egtak_RB_std = std(interclusterdist_egtak_RB);
interclusterdist_egtak_GB_std = std(interclusterdist_egtak_GB);


clustersize_RC_mean = [size_phys_RC_mean size_egta_RC_mean size_egtak_RC_mean];
clustersize_GC_mean = [size_phys_GC_mean size_egta_GC_mean size_egtak_GC_mean];
clustersize_BC_mean = [size_phys_BC_mean size_egta_BC_mean size_egtak_BC_mean];
interclusterdist_RG_mean = [interclusterdist_phys_RG_mean interclusterdist_egta_RG_mean interclusterdist_egtak_RG_mean];
interclusterdist_RB_mean = [interclusterdist_phys_RB_mean interclusterdist_egta_RB_mean interclusterdist_egtak_RB_mean];
interclusterdist_GB_mean = [interclusterdist_phys_GB_mean interclusterdist_egta_GB_mean interclusterdist_egtak_GB_mean];

clustersize_RC_std = [size_phys_RC_std size_egta_RC_std size_egtak_RC_std];
clustersize_GC_std = [size_phys_GC_std size_egta_GC_std size_egtak_GC_std];
clustersize_BC_std = [size_phys_BC_std size_egta_BC_std size_egtak_BC_std];
interclusterdist_RG_std = [interclusterdist_phys_RG_std interclusterdist_egta_RG_std interclusterdist_egtak_RG_std];
interclusterdist_RB_std = [interclusterdist_phys_RB_std interclusterdist_egta_RB_std interclusterdist_egtak_RB_std];
interclusterdist_GB_std = [interclusterdist_phys_GB_std interclusterdist_egta_GB_std interclusterdist_egtak_GB_std];


% Rescale the data to fit on a plot (too many differences in scales otherwise)
norm_area_mCling_mean              = area_mCling_mean/area_mCling_mean(1);
norm_overlapWithGreen_mean         = overlapWithGreen_mean/overlapWithGreen_mean(1);
norm_overlapWithBlue_mean          = overlapWithBlue_mean/overlapWithBlue_mean(1);
norm_weightedOverlapWithGreen_mean = weightedOverlapWithGreen_mean/weightedOverlapWithGreen_mean(1);
norm_weightedOverlapWithBlue_mean  = weightedOverlapWithBlue_mean/weightedOverlapWithBlue_mean(1);
norm_clustersize_RC_mean           = clustersize_RC_mean/clustersize_RC_mean(1);
norm_clustersize_GC_mean           = clustersize_GC_mean/clustersize_GC_mean(1);
norm_clustersize_BC_mean           = clustersize_BC_mean/clustersize_BC_mean(1);
norm_interclusterdist_RG_mean      = interclusterdist_RG_mean/interclusterdist_RG_mean(1);
norm_interclusterdist_RB_mean      = interclusterdist_RB_mean/interclusterdist_RB_mean(1);
norm_interclusterdist_GB_mean      = interclusterdist_GB_mean/interclusterdist_GB_mean(1);

norm_area_mCling_std               = area_mCling_std/area_mCling_mean(1);
norm_overlapWithGreen_std          = overlapWithGreen_std/overlapWithGreen_mean(1);
norm_overlapWithBlue_std           = overlapWithBlue_std/overlapWithBlue_mean(1);
norm_weightedOverlapWithGreen_std  = weightedOverlapWithGreen_std/weightedOverlapWithGreen_mean(1);
norm_weightedOverlapWithBlue_std   = weightedOverlapWithBlue_std/weightedOverlapWithBlue_mean(1);
norm_clustersize_RC_std            = clustersize_RC_std/clustersize_RC_mean(1);
norm_clustersize_GC_std            = clustersize_GC_std/clustersize_GC_mean(1);
norm_clustersize_BC_std            = clustersize_BC_std/clustersize_BC_mean(1);
norm_interclusterdist_RG_std       = interclusterdist_RG_std/interclusterdist_RG_mean(1);
norm_interclusterdist_RB_std       = interclusterdist_RB_std/interclusterdist_RB_mean(1);
norm_interclusterdist_GB_std       = interclusterdist_GB_std/interclusterdist_GB_mean(1);


fig1 = figure;
clear var ctr ydt
data = [norm_area_mCling_mean' norm_overlapWithGreen_mean' ...
        norm_overlapWithBlue_mean' norm_weightedOverlapWithGreen_mean' ...
        norm_weightedOverlapWithBlue_mean'];
Std  = [norm_area_mCling_std' norm_overlapWithGreen_std' ...
        norm_overlapWithBlue_std' norm_weightedOverlapWithGreen_std' ...
        norm_weightedOverlapWithBlue_std'];
hBar = bar(1:3,data);
for i = 1:size(data,2)
    ctr(i,:) = bsxfun(@plus, hBar(1).XData, [hBar(i).XOffset]');
    ydt(i,:) = hBar(i).YData;
end
hold on
errorbar(ctr, ydt, Std', '.k')
set(gca,'XTick',1:3,'XTickLabel',{'PHYS','EGTA','EGTA/K+'})
legend({'Area mCling','Overlap mCling/a-syn','Overlap mCling/VAMP2','Weighted overlap mCling/a-syn','Weighted overlap mCling/VAMP2'},'Location','northwest')
ylabel('Value normalized by value PHYS condition')
title(strcat('Influence of area on overlap (',temperature,' °C)'));
set(gca,'fontsize',14);
hold off


fig2 = figure;
clear var ctr ydt
data = [norm_area_mCling_mean' norm_clustersize_RC_mean' ...
        norm_clustersize_GC_mean' norm_clustersize_BC_mean' ...
        norm_overlapWithGreen_mean' norm_overlapWithBlue_mean' ...
        norm_weightedOverlapWithGreen_mean' norm_weightedOverlapWithBlue_mean'];
Std  = [norm_area_mCling_std' norm_clustersize_RC_std' ...
        norm_clustersize_GC_std' norm_clustersize_BC_std' ...
        norm_overlapWithGreen_std' norm_overlapWithBlue_std' ...
        norm_weightedOverlapWithGreen_std' norm_weightedOverlapWithBlue_std'];
hBar = bar(1:3,data);
for i = 1:size(data,2)
    ctr(i,:) = bsxfun(@plus, hBar(1).XData, [hBar(i).XOffset]');
    ydt(i,:) = hBar(i).YData;
end
hold on
errorbar(ctr, ydt, Std', '.k')
set(gca,'XTick',1:3,'XTickLabel',{'PHYS','EGTA','EGTA/K+'})
legend({'Area mCling','Size mCling','Size a-synuclein','Size VAMP2','Overlap mCling/a-syn','Overlap mCling/VAMP2','Weighted overlap mCling/a-syn','Weighted overlap mCling/VAMP2'},'Location','northwest')
ylabel('Value normalized by value PHYS condition')
title(strcat('Influence of area on overlap (',temperature,' °C)'));
set(gca,'fontsize',14);


fig3 = figure;
clear var ctr ydt
data = [interclusterdist_RG_mean' interclusterdist_RB_mean' interclusterdist_GB_mean'];
Std  = [interclusterdist_RG_std'  interclusterdist_RB_std'  interclusterdist_GB_std'];
hBar = bar(1:3,data);
for i = 1:size(data,2)
    ctr(i,:) = bsxfun(@plus, hBar(1).XData, [hBar(i).XOffset]');
    ydt(i,:) = hBar(i).YData;
end
hold on
errorbar(ctr, ydt, Std', '.k')
set(gca,'XTick',1:3,'XTickLabel',{'PHYS','EGTA','EGTA/K+'})
legend({'Intercluster distance mCling/a-syn','Intercluster distance mCling/VAMP2','Intercluster distance a-syn/VAMP2'},'Location','northwest')
ylabel('Value normalized by value PHYS condition')
title(strcat('Intercluster distance comparison (',temperature,' °C)'));
set(gca,'fontsize',14);
hold off

fig4 = figure;
clear var ctr ydt
data = [norm_weightedOverlapWithGreen_mean' norm_weightedOverlapWithBlue_mean'];
Std  = [norm_weightedOverlapWithGreen_std'  norm_weightedOverlapWithBlue_std'];
hBar = bar(1:3,data);
for i = 1:size(data,2)
    ctr(i,:) = bsxfun(@plus, hBar(1).XData, [hBar(i).XOffset]');
    ydt(i,:) = hBar(i).YData;
end
hold on
errorbar(ctr, ydt, Std', '.k')
set(gca,'XTick',1:3,'XTickLabel',{'PHYS','EGTA','EGTA/K+'})
legend({'Weighted overlap mCling/a-syn','Weighted overlap mCling/VAMP2'},'Location','northwest')
ylabel('Value normalized by value PHYS condition')
title(strcat('Weighted overlap comparison (',temperature,' °C)'));
set(gca,'fontsize',14);
hold off

path_fig1 = fullfile(path_output,strcat('groupedBarplot_overlap_',temperature,'C-1.fig'));
path_fig2 = fullfile(path_output,strcat('groupedBarplot_overlap_',temperature,'C-2.fig'));
path_fig3 = fullfile(path_output,strcat('groupedBarplot_interactiondist_',temperature,'C.fig'));
path_fig4 = fullfile(path_output,strcat('groupedBarplot_weighted_overlaps_',temperature,'C.fig'));

savefig(fig1,path_fig1)
savefig(fig2,path_fig2)
savefig(fig3,path_fig3)
savefig(fig4,path_fig4)
