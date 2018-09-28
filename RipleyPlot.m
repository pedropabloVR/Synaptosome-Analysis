function fig = RipleyPlot(r_hist,H_all,color,title_plot,ylimits)

if strcmp(color,'red')
    light = [1 0.8 0.8];
    dark = [0.7 0 0];
elseif strcmp(color,'green')
    light = [0.8 1 0.8];
    dark = [0 0.7 0];
elseif strcmp(color,'blue')
    light = [0.8 0.8 1];
    dark = [0 0 0.7];   
else
    disp('Pick a valid colour for RipleyPlot.')
end
    
fig = figure;

% Ripley curves for all the ROIs
h1 = plot(r_hist,H_all,'Color',light);
hold on

% Mean of the Ripley curves of all the ROIs
H_all_mean = mean(H_all,2);
h2 = plot(r_hist,H_all_mean,'Color',dark,'LineWidth',1.5);
[idx_x,~] = getCoordinatesMax(r_hist,H_all_mean);
l = line([idx_x idx_x], [0 max(H_all_mean)],'Color','black','LineStyle','--','LineWidth',1.5);
hold on

xlabel('r (nm)');
ylabel('Ripley''s L(r) - r function');
title(title_plot);
legend([h1(1) h2 l],{'Data','Mean','Estimate radius'},'Location','southeast');
set(gca,'fontsize',14);
ylim(ylimits)

end