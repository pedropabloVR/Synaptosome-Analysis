% This function gets the means and standard deviations for data that has
% been seperated into 3 conditions (specified by condition_levels in a
% second column).
% 
% INPUT:
%   X...............: column with data
%   condition_col...: column with labels
%   condition_levels: values of the labels (e.g. {'phys','egta','egtak'})
% 
% OUTPUT:
%   Mean............: 1x3 array with means of the three conditions
%   STD.............: 1x3 array with sd's of the three conditions
% 
% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated on 30 Sept 2018


function [Mean,STD] = getMeanAndStdConditions(X,condition_col,condition_levels)

% Get names of the conditions
cond1 = condition_levels{1};
cond2 = condition_levels{2};
cond3 = condition_levels{3};

% Calculate the mean for each condition
X_mean_cond1 = mean(X(strcmp(condition_col,cond1)));
X_mean_cond2 = mean(X(strcmp(condition_col,cond2)));
X_mean_cond3 = mean(X(strcmp(condition_col,cond3)));

% Calculate the standard deviation for each condition
X_std_cond1  = std(X(strcmp(condition_col,cond1)));
X_std_cond2  = std(X(strcmp(condition_col,cond2)));
X_std_cond3  = std(X(strcmp(condition_col,cond3)));

% Put the means and standard deviations of the three conditions together
Mean = [X_mean_cond1 X_mean_cond2 X_mean_cond3];
STD = [X_std_cond1 X_std_cond2 X_std_cond3];

end