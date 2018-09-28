function [Mean,STD] = getMeanAndStdConditions(X,condition_col,condition_levels)

cond1 = condition_levels{1};
cond2 = condition_levels{2};
cond3 = condition_levels{3};

X_mean_cond1 = mean(X(strcmp(condition_col,cond1)));
X_mean_cond2 = mean(X(strcmp(condition_col,cond2)));
X_mean_cond3 = mean(X(strcmp(condition_col,cond3)));

X_std_cond1  = std(X(strcmp(condition_col,cond1)));
X_std_cond2  = std(X(strcmp(condition_col,cond2)));
X_std_cond3  = std(X(strcmp(condition_col,cond3)));

Mean = [X_mean_cond1 X_mean_cond2 X_mean_cond3];
STD = [X_std_cond1 X_std_cond2 X_std_cond3];

end