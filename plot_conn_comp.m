num_components_smooth = movmean(num_components,100);
mean_comp_cm = (component_mean .* 400) ./ 100000000;
mean_comp_smooth = movmean(mean_comp_cm, 100);
figure();
yyaxis left
plot(depth(939:3392),num_components(939:3392), 'Color',[0 0.5 1]);
hold on
plot(depth(939:3392),num_components_smooth(939:3392),'LineWidth',5,'LineStyle','-','Color',[0 0.4470 0.7410]);
xlabel('Distance from top (cm)');
xlim([1.5 7])
ylabel('Number of components');
yyaxis right
plot(depth(939:3392), mean_comp_cm(939:3392), 'Color',[1 0.5 0]);
hold on
plot(depth(939:3392), mean_comp_smooth(939:3392),'LineWidth',5,'LineStyle','-','Color',[0.8500 0.3250 0.0980]);
ylabel('Mean component area (cm^2)');
title('Connected Components of SM 117.71');
corr_coef = corrcoef(num_components_smooth(939:3392),mean_comp_smooth(939:3392));
disp(corr_coef(2,1));