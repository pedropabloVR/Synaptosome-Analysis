function fig = RipleyPlotMeansConditions(r_hist,H_all_1,H_all_2,H_all_3,plot_legend,title_plot)

fig = figure;

% Condition 1
H_all_mean1 = mean(H_all_1,2);
h1 = plot(r_hist,H_all_mean1,'-k','LineWidth',1.5);
[idx_x,~] = getCoordinatesMax(r_hist,H_all_mean1);
line([idx_x idx_x], [0 max(H_all_mean1)],'Color','black','LineStyle','-','LineWidth',1.5);
hold on

% Condition 2
H_all_mean2 = mean(H_all_2,2);
h2 = plot(r_hist,H_all_mean2,'--k','LineWidth',1.5);
[idx_x,~] = getCoordinatesMax(r_hist,H_all_mean2);
line([idx_x idx_x], [0 max(H_all_mean2)],'Color','black','LineStyle','--','LineWidth',1.5);
hold on

% Condition 3
H_all_mean3 = mean(H_all_3,2);
h3 = plot(r_hist,H_all_mean3,':','LineWidth',1.5);
[idx_x,~] = getCoordinatesMax(r_hist,H_all_mean3);
line([idx_x idx_x], [0 max(H_all_mean3)],'Color','black','LineStyle',':','LineWidth',1.5);
hold off

legend([h1,h2,h3],plot_legend,'Location','southeast')

xlabel('r (nm)');
ylabel('Ripley''s L(r) - r function');
title(title_plot);
set(gca,'fontsize',14);
ylim([0 1.1*max([H_all_mean1; H_all_mean2; H_all_mean3])])

end