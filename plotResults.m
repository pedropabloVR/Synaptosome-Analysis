function [] = plotResults(path_output, results_PHYS, results_EGTA, ...
    results_EGTAK, magnification, flagprint)

%% Plot and compare results from overlap analysis

path_figures = fullfile(path_output,'figures');
mkdir(path_figures);
path_figures_overlap = fullfile(path_figures,'overlap_analysis');
mkdir(path_figures_overlap);

% Groups
[condition_PHYS{1:size(results_PHYS,1)}]   = deal('PHYS');
[condition_EGTA{1:size(results_EGTA,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(results_EGTAK,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];


% Area of synaptosomes ----------------------------------------------------

areaRed_PHYS = results_PHYS.Area;
areaRed_EGTA = results_EGTA.Area;
areaRed_EGTAK = results_EGTAK.Area;
areaRed = [areaRed_PHYS; areaRed_EGTA; areaRed_EGTAK];

% Beeswarm and boxplot
fig1 = figure;
subplot(121)
h = plotSpread({areaRed_PHYS*magnification^2,areaRed_EGTA*magnification^2,areaRed_EGTAK*magnification^2},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Area (nm^2)');
title('Synaptosome area');
pbaspect([1 1 1])

subplot(122)
boxplot(areaRed*magnification^2,conditions);
ylabel('Area (nm^2)');
title('Synaptosome area');
pbaspect([1 1 1])

if flagprint
    savefig(fig1,fullfile(path_figures_overlap,'area_synaptosomes.fig'))
end


% Area of a-synuclein -----------------------------------------------------

areaGreen_PHYS = results_PHYS.AreaGreen;
areaGreen_EGTA = results_EGTA.AreaGreen;
areaGreen_EGTAK = results_EGTAK.AreaGreen;
areaGreen = [areaGreen_PHYS; areaGreen_EGTA; areaGreen_EGTAK];

% Beeswarm and boxplot
fig4 = figure;
subplot(121)
h = plotSpread({areaGreen_PHYS*magnification^2,areaGreen_EGTA*magnification^2,areaGreen_EGTAK*magnification^2},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Area (nm^2)');
title('a-synuclein area');
pbaspect([1 1 1])

subplot(122)
boxplot(areaGreen*magnification^2,conditions);
ylabel('Area (nm^2)');
title('a-synuclein area');
pbaspect([1 1 1])

if flagprint
    savefig(fig4,fullfile(path_figures_overlap,'area-a-synuclein.fig'))
end


% Area of VAMP2 -----------------------------------------------------------

areaBlue_PHYS = results_PHYS.AreaBlue;
areaBlue_EGTA = results_EGTA.AreaBlue;
areaBlue_EGTAK = results_EGTAK.AreaBlue;
areaBlue = [areaBlue_PHYS; areaBlue_EGTA; areaBlue_EGTAK];

% Beeswarm and boxplot
fig5 = figure;
subplot(121)
h = plotSpread({areaBlue_PHYS*magnification^2,areaBlue_EGTA*magnification^2,areaBlue_EGTAK*magnification^2},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Area (nm^2)');
title('VAMP2 area');
pbaspect([1 1 1])

subplot(122)
boxplot(areaBlue*magnification^2,conditions);
ylabel('Area (nm^2)');
title('VAMP2 area');
pbaspect([1 1 1])

if flagprint
    savefig(fig5,fullfile(path_figures_overlap,'area-VAMP2.fig'))
end


% Overlap of synaptosomes with a-synuclein --------------------------------

overlapGreen_PHYS = results_PHYS.OverlapWithGreen;
overlapGreen_EGTA = results_EGTA.OverlapWithGreen;
overlapGreen_EGTAK = results_EGTAK.OverlapWithGreen;
overlapGreen = [overlapGreen_PHYS; overlapGreen_EGTA; overlapGreen_EGTAK];

% Beeswarm and boxplot
fig2 = figure;
subplot(121)
h = plotSpread({overlapGreen_PHYS,overlapGreen_EGTA,overlapGreen_EGTAK},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('% area overlap'); ylim([0 100]);
title({'Unweighted overlap','mCling and a-synuclein'});
pbaspect([1 1 1])
ylim([0 100])

subplot(122)
boxplot(overlapGreen,conditions);
ylabel('% area overlap'); ylim([0 100]);
title({'Unweighted overlap','mCling and a-synuclein'});
pbaspect([1 1 1])
ylim([0 100])

if flagprint
    savefig(fig2,fullfile(path_figures_overlap,'unweighted-overlap-mcling-a-synuclein.fig'))
    % export figure to a pdf for display
    filename = fullfile(path_figures_overlap,'unweighted-overlap-mcling-a-synuclein');
    print(fig2,filename,'-dpdf','-r300')
end


% Overlap of synaptosomes with VAMP2 --------------------------------------

overlapBlue_PHYS = results_PHYS.OverlapWithBlue;
overlapBlue_EGTA = results_EGTA.OverlapWithBlue;
overlapBlue_EGTAK = results_EGTAK.OverlapWithBlue;
overlapBlue = [overlapBlue_PHYS; overlapBlue_EGTA; overlapBlue_EGTAK];

% Beeswarm and boxplot
fig3 = figure;
subplot(121)
h = plotSpread({overlapBlue_PHYS,overlapBlue_EGTA,overlapBlue_EGTAK},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('% area overlap'); ylim([0 100]);
title({'Unweighted overlap','mCling and VAMP2'});
pbaspect([1 1 1])
ylim([0 100])

subplot(122)
boxplot(overlapBlue,conditions);
ylabel('% area overlap'); ylim([0 100]);
title({'Unweighted overlap','mCling and VAMP2'});
pbaspect([1 1 1])
ylim([0 100])

if flagprint
    savefig(fig3,fullfile(path_figures_overlap,'unweighted-overlap-mcling-vamp2.fig'))
    filename = fullfile(path_figures_overlap,'unweighted-overlap-mcling-vamp2');
    print(fig3,filename,'-dpdf','-r300')
end


% Weighted overlap of synaptosomes with a-synuclein -----------------------

weightedOverlapWithGreen_PHYS  = results_PHYS.WeightedOverlapWithGreen;
weightedOverlapWithGreen_EGTA  = results_EGTA.WeightedOverlapWithGreen;
weightedOverlapWithGreen_EGTAK = results_EGTAK.WeightedOverlapWithGreen;
weightedOverlapWithGreen = [weightedOverlapWithGreen_PHYS;
                            weightedOverlapWithGreen_EGTA;
                            weightedOverlapWithGreen_EGTAK];
fig6 = figure;
subplot(121)
h = plotSpread({weightedOverlapWithGreen_PHYS,weightedOverlapWithGreen_EGTA,weightedOverlapWithGreen_EGTAK},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Weighted overlap');
title({'Weighted overlap','mCling and a-synuclein'});
pbaspect([1 1 1])
ylim([0 1])

subplot(122)
boxplot(weightedOverlapWithGreen,conditions);
ylabel('Weighted overlap');
title({'Weighted overlap','mCling and a-synuclein'});
pbaspect([1 1 1])
ylim([0 1])

if flagprint
    savefig(fig6,fullfile(path_figures_overlap,'weighted-overlap-mcling-a-synuclein.fig'))
    filename = fullfile(path_figures_overlap,'weighted-overlap-mcling-a-synuclein');
    print(fig6,filename,'-dpdf','-r300')
end


% Weighted overlap of synaptosomes with VAMP2 -----------------------------

weightedOverlapWithBlue_PHYS  = results_PHYS.WeightedOverlapWithBlue;
weightedOverlapWithBlue_EGTA  = results_EGTA.WeightedOverlapWithBlue;
weightedOverlapWithBlue_EGTAK = results_EGTAK.WeightedOverlapWithBlue;
weightedOverlapWithBlue = [weightedOverlapWithBlue_PHYS;
                           weightedOverlapWithBlue_EGTA;
                           weightedOverlapWithBlue_EGTAK];
fig7 = figure;
subplot(121)
h = plotSpread({weightedOverlapWithBlue_PHYS,weightedOverlapWithBlue_EGTA,weightedOverlapWithBlue_EGTAK},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Weighted overlap');
title({'Weighted overlap','mCling and VAMP2'});
pbaspect([1 1 1])
ylim([0 1])

subplot(122)
boxplot(weightedOverlapWithBlue,conditions);
ylabel('Weighted overlap');
title({'Weighted overlap','mCling and VAMP2'});
pbaspect([1 1 1])
ylim([0 1])

if flagprint
    savefig(fig7,fullfile(path_figures_overlap,'weighted-overlap-mcling-vamp2.fig'))
    filename = fullfile(path_figures_overlap,'weighted-overlap-mcling-vamp2');
    print(fig7,filename,'-dpdf','-r300')
end

%% Plot and compare results from Ripley analysis

path_figures_ripley = fullfile(path_figures,'ripley_analysis');
mkdir(path_figures_ripley);

% Clustersize mCling from Ripley ------------------------------------------

clustersize_phys_RC = results_PHYS.clustersizeRC;
clustersize_egta_RC = results_EGTA.clustersizeRC;
clustersize_egtak_RC = results_EGTAK.clustersizeRC;

% Get groups
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(clustersize_phys_RC,1)}]   = deal('PHYS');
[condition_EGTA{1:size(clustersize_egta_RC,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(clustersize_egtak_RC,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot clustersizes estimated from Ripley's curves
fig8 = figure;

subplot(121) % beeswarm plot
h = plotSpread({clustersize_phys_RC,clustersize_egta_RC,clustersize_egtak_RC},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Cluster size (nm)');
title('Clustersize mCling');
set(gca,'fontsize',14);
%ylim([0 1000])
pbaspect([1 1 1])

subplot(122) % boxplot
boxplot([clustersize_phys_RC; clustersize_egta_RC; clustersize_egtak_RC],conditions);
ylabel('Cluster size (nm)');
title('Clustersize mCling');
set(gca,'fontsize',14);
%ylim([0 1000])
pbaspect([1 1 1])

if flagprint
    savefig(fig8,fullfile(path_figures_ripley,'clustersize_mCling.fig'))
end

% One-way ANOVA
p = anova1([clustersize_phys_RC; clustersize_egta_RC; clustersize_egtak_RC],conditions,'off');
if p < 0.05
    disp(['Significant difference in mCling clustersize between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in mCling clustersize between conditions (p-value = ' num2str(p) ')'])
end

% Clustersize a-synuclein from Ripley -------------------------------------

clustersize_phys_GC = results_PHYS.clustersizeGC;
clustersize_egta_GC = results_EGTA.clustersizeGC;
clustersize_egtak_GC = results_EGTAK.clustersizeGC;

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(clustersize_phys_GC,1)}]   = deal('PHYS');
[condition_EGTA{1:size(clustersize_egta_GC,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(clustersize_egtak_GC,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot clustersizes estimated from Ripley's curves
fig9 = figure;
subplot(121) % beeswarm plot
h = plotSpread({clustersize_phys_GC,clustersize_egta_GC,clustersize_egtak_GC},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Cluster size (nm)');
title('Clustersize a-synuclein');
set(gca,'fontsize',14);
%ylim([0 1000])
pbaspect([1 1 1])

subplot(122) % boxplot
boxplot([clustersize_phys_GC; clustersize_egta_GC; clustersize_egtak_GC],conditions);
ylabel('Cluster size (nm)');
title('Clustersize a-synuclein');
set(gca,'fontsize',14);
%ylim([0 1000])
pbaspect([1 1 1])

if flagprint
    savefig(fig9,fullfile(path_figures_ripley,'clustersize_a-synuclein.fig'))
end

% One-way ANOVA
p = anova1([clustersize_phys_GC; clustersize_egta_GC; clustersize_egtak_GC],conditions,'off');
if p < 0.05
    disp(['Significant difference in a-synuclein clustersize between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in a-synuclein clustersize between conditions (p-value = ' num2str(p) ')'])
end

% Clustersize VAMP2 from Ripley -------------------------------------------

clustersize_phys_BC = results_PHYS.clustersizeBC;
clustersize_egta_BC = results_EGTA.clustersizeBC;
clustersize_egtak_BC = results_EGTAK.clustersizeBC;

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(clustersize_phys_BC,1)}]   = deal('PHYS');
[condition_EGTA{1:size(clustersize_egta_BC,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(clustersize_egtak_BC,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot clustersizes estimated from Ripley's curves
fig10 = figure;
subplot(121) % beeswarm plot
h = plotSpread({clustersize_phys_BC,clustersize_egta_BC,clustersize_egtak_BC},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Cluster size (nm)');
title('Clustersize VAMP2');
set(gca,'fontsize',14);
%ylim([0 1000])
pbaspect([1 1 1])

subplot(122) % boxplot
boxplot([clustersize_phys_BC; clustersize_egta_BC; clustersize_egtak_BC],conditions);
ylabel('Cluster size (nm)');
title('Clustersize VAMP2');
set(gca,'fontsize',14);
%ylim([0 1000])
pbaspect([1 1 1])

if flagprint
    savefig(fig10,fullfile(path_figures_ripley,'clustersize_VAMP2.fig'))
end

% One-way ANOVA
p = anova1([clustersize_phys_BC; clustersize_egta_BC; clustersize_egtak_BC],conditions,'off');
if p < 0.05
    disp(['Significant difference in VAMP2 clustersize between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in VAMP2 clustersize between conditions (p-value = ' num2str(p) ')'])
end


% Intercluster distance between mCling and a-synuclein from Ripley --------

interclusterdist_phys_RG = results_PHYS.interclusterdistRG;
interclusterdist_egta_RG = results_EGTA.interclusterdistRG;
interclusterdist_egtak_RG = results_EGTAK.interclusterdistRG;

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(interclusterdist_phys_RG,1)}]   = deal('PHYS');
[condition_EGTA{1:size(interclusterdist_egta_RG,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(interclusterdist_egtak_RG,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot intercluster distances estimated from Ripley's curves
fig11 = figure;
subplot(121) % beeswarm plot
h = plotSpread({interclusterdist_phys_RG,interclusterdist_egta_RG,interclusterdist_egtak_RG},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','mCling and a-synuclein'});
set(gca,'fontsize',14);
pbaspect([1 1 1])

subplot(122) % boxplot
boxplot([interclusterdist_phys_RG; interclusterdist_egta_RG; interclusterdist_egtak_RG],conditions);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','mCling and a-synuclein'});
set(gca,'fontsize',14);
pbaspect([1 1 1])

if flagprint
    savefig(fig11,fullfile(path_figures_ripley,'intercluster_distance_mCling_a-synuclein.fig'))
end

% One-way ANOVA
p = anova1([interclusterdist_phys_RG; interclusterdist_egta_RG; interclusterdist_egtak_RG],conditions,'off');
if p < 0.05
    disp(['Significant difference in mCling/a-synuclein intercluster distance between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in mCling/a-synuclein intercluster distance between conditions (p-value = ' num2str(p) ')'])
end


% Intercluster distance between mCling and VAMP2 from Ripley --------------

interclusterdist_phys_RB = results_PHYS.interclusterdistRB;
interclusterdist_egta_RB = results_EGTA.interclusterdistRB;
interclusterdist_egtak_RB = results_EGTAK.interclusterdistRB;

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(interclusterdist_phys_RB,1)}]   = deal('PHYS');
[condition_EGTA{1:size(interclusterdist_egta_RB,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(interclusterdist_egtak_RB,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot intercluster distances estimated from Ripley's curves
fig12 = figure;
subplot(121) % beeswarm plot
h = plotSpread({interclusterdist_phys_RB,interclusterdist_egta_RB,interclusterdist_egtak_RB},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','mCling and VAMP2'});
set(gca,'fontsize',14);
pbaspect([1 1 1])

subplot(122) % boxplot
boxplot([interclusterdist_phys_RB; interclusterdist_egta_RB; interclusterdist_egtak_RB],conditions);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','mCling and VAMP2'});
set(gca,'fontsize',14);
pbaspect([1 1 1])

if flagprint
    savefig(fig12,fullfile(path_figures_ripley,'intercluster_distance_mCling_VAMP2.fig'))
end

% One-way ANOVA
p = anova1([interclusterdist_phys_RB; interclusterdist_egta_RB; interclusterdist_egtak_RB],conditions,'off');
if p < 0.05
    disp(['Significant difference in mCling/VAMP2 intercluster distance between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in mCling/VAMP2 intercluster distance between conditions (p-value = ' num2str(p) ')'])
end


% Intercluster distance between a-synuclein and VAMP2 from Ripley ---------

interclusterdist_phys_GB = results_PHYS.interclusterdistGB;
interclusterdist_egta_GB = results_EGTA.interclusterdistGB;
interclusterdist_egtak_GB = results_EGTAK.interclusterdistGB;

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(interclusterdist_phys_GB,1)}]   = deal('PHYS');
[condition_EGTA{1:size(interclusterdist_egta_GB,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(interclusterdist_egtak_GB,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot intercluster distances estimated from Ripley's curves
fig13 = figure;
subplot(121) % beeswarm plot
h = plotSpread({interclusterdist_phys_GB,interclusterdist_egta_GB,interclusterdist_egtak_GB},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','a-synuclein and VAMP2'});
set(gca,'fontsize',14);
pbaspect([1 1 1])

subplot(122) % boxplot
boxplot([interclusterdist_phys_GB; interclusterdist_egta_GB; interclusterdist_egtak_GB],conditions);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','a-synuclein and VAMP2'});
set(gca,'fontsize',14);
pbaspect([1 1 1])

if flagprint
    savefig(fig13,fullfile(path_figures_ripley,'intercluster_distance_a-synuclein_VAMP2.fig'))
end

% One-way ANOVA
p = anova1([interclusterdist_phys_GB; interclusterdist_egta_GB; interclusterdist_egtak_GB],conditions,'off');
if p < 0.05
    disp(['Significant difference in a-synuclein/VAMP2 intercluster distance between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in a-synuclein/VAMP2 intercluster distance between conditions (p-value = ' num2str(p) ')'])
end

% %% Plot the Ripley's functions
% 
% path_ripley_curves = fullfile(path_figures_ripley,'curves');
% mkdir(path_ripley_curves);
% 
% ymin = 0;
% ymax = 6000;
% 
% % mCLING
% 
% 
% fig1 = RipleyPlot(r_hist,H_all_phys_RC(:,all(~isnan(H_all_phys_RC))),'red','mCling PHYS',[ymin ymax]);
% fig2 = RipleyPlot(r_hist,H_all_egta_RC(:,all(~isnan(H_all_egta_RC))),'red','mCling EGTA',[ymin ymax]);
% fig3 = RipleyPlot(r_hist,H_all_egtak_RC(:,all(~isnan(H_all_egtak_RC))),'red','mCling EGTA/K+',[ymin ymax]);
% fig4 = RipleyPlotMeansConditions(r_hist, H_all_phys_RC(:,all(~isnan(H_all_phys_RC))),...
%                                          H_all_egta_RC(:,all(~isnan(H_all_egta_RC))),...
%                                          H_all_egtak_RC(:,all(~isnan(H_all_egtak_RC))),...
%                                          {'PHYS','EGTA','EGTA/K+'},'Mean of Ripley''s - mCling');
% 
% % a-synuclein
% fig5 = RipleyPlot(r_hist,H_all_phys_GC(:,all(~isnan(H_all_phys_GC))), 'green','a-synuclein PHYS',[ymin ymax]);
% fig6 = RipleyPlot(r_hist,H_all_egta_GC(:,all(~isnan(H_all_egta_GC))), 'green','a-synuclein EGTA',[ymin ymax]);
% fig7 = RipleyPlot(r_hist,H_all_egtak_GC(:,all(~isnan(H_all_egtak_GC))),'green','a-synuclein EGTA/K+',[ymin ymax]);
% fig8 = RipleyPlotMeansConditions(r_hist, H_all_phys_GC(:,all(~isnan(H_all_phys_GC))),...
%                                          H_all_egta_GC(:,all(~isnan(H_all_egta_GC))),...
%                                          H_all_egtak_GC(:,all(~isnan(H_all_egtak_GC))),...
%                                          {'PHYS','EGTA','EGTA/K+'},'Mean of Ripley''s - a-synuclein');
% 
% % VAMP2
% fig9  = RipleyPlot(r_hist,H_all_phys_BC(:,all(~isnan(H_all_phys_BC))), 'blue','VAMP2 PHYS',[ymin ymax]);
% fig10 = RipleyPlot(r_hist,H_all_egta_BC(:,all(~isnan(H_all_egta_BC))), 'blue','VAMP2 EGTA',[ymin ymax]);
% fig11 = RipleyPlot(r_hist,H_all_egtak_BC(:,all(~isnan(H_all_egtak_BC))),'blue','VAMP2 EGTA/K+',[ymin ymax]);
% fig12 = RipleyPlotMeansConditions(r_hist, H_all_phys_BC(:,all(~isnan(H_all_phys_BC))),...
%                                           H_all_egta_BC(:,all(~isnan(H_all_egta_BC))),...
%                                           H_all_egtak_BC(:,all(~isnan(H_all_egtak_BC))),...
%                                           {'PHYS','EGTA','EGTA/K+'},'Mean of Ripley''s - VAMP2');
% 
% if flagprint
%     savefig(fig1 ,fullfile(path_ripley_curves,'ripley_mcling_phys.fig'))
%     savefig(fig2 ,fullfile(path_ripley_curves,'ripley_mcling_egta.fig'))
%     savefig(fig3 ,fullfile(path_ripley_curves,'ripley_mcling_egtak.fig'))
%     savefig(fig4 ,fullfile(path_ripley_curves,'ripley_mcling.fig'))
%     savefig(fig5 ,fullfile(path_ripley_curves,'ripley_a-synuclein_phys.fig'))
%     savefig(fig6 ,fullfile(path_ripley_curves,'ripley_a-synuclein_egta.fig'))
%     savefig(fig7 ,fullfile(path_ripley_curves,'ripley_a-synuclein_egtak.fig'))
%     savefig(fig8 ,fullfile(path_ripley_curves,'ripley_a-synuclein.fig'))
%     savefig(fig9 ,fullfile(path_ripley_curves,'ripley_VAMP2_phys.fig'))
%     savefig(fig10,fullfile(path_ripley_curves,'ripley_VAMP2_egta.fig'))
%     savefig(fig11,fullfile(path_ripley_curves,'ripley_VAMP2_egtak.fig'))
%     savefig(fig12,fullfile(path_ripley_curves,'ripley_VAMP2.fig'))
% end

end
