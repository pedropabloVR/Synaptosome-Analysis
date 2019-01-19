
% Script to plot the comparison between synaptosome conditions using the
% output file from script compare_conditions.m
% 15/01/2019

% Import .csv file with results_combined_pooled after manual deletion of
% outliers

repeats         = {'A','B'}; %,'Oct-Nov'};
output_dir      = fullfile('E:\Experiments\synaptosomes\analysis_20190107\37C_results_thresh30\Results_combined');
magnification   = 10;
flagprint       = 1; % set to 1 to save the fig visualization of the results

%% Write away pooled results from all the repeats in a new folder
results_pooled = resultscombinedpooledtrimmed;

% Create new output folder
path_output = fullfile(output_dir,'Results_combined_pooled_culled');
mkdir(path_output);

% Split results table into conditions (for convenience when plotting)
results_PHYS = results_pooled(strcmp(results_pooled.condition,'phys'),:);
results_EGTA = results_pooled(strcmp(results_pooled.condition,'egta'),:);
results_EGTAK = results_pooled(strcmp(results_pooled.condition,'egtak'),:);

% Plot all results
plotResults(path_output,results_PHYS,results_EGTA,results_EGTAK,magnification,flagprint)


