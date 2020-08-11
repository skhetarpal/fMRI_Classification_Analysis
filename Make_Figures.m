for i=1:length(x)
    age(i) = x(i).age/30;
    %Collect data for brain region V1
    acc1(i) = x(i).apAnalyses(3).apROI{1}.acc;
    pval1(i) = x(i).apAnalyses(3).apROI{1}.pval;
    perm_upper1(i) = x(i).apAnalyses(3).apROI{1}.perm_upper;
    num_obs1(i) = x(i).apAnalyses(3).apROI{1}.num_obs;
    num_vox1(i) = x(i).apAnalyses(3).apROI{1}.num_vox;
    %Collect data for brain region PIT    
    acc2(i) = x(i).apAnalyses(3).apROI{7}.acc;
    pval2(i) = x(i).apAnalyses(3).apROI{7}.pval;
    perm_upper2(i) = x(i).apAnalyses(3).apROI{7}.perm_upper;
    num_obs2(i) = x(i).apAnalyses(3).apROI{7}.num_obs;
    num_vox2(i) = x(i).apAnalyses(3).apROI{7}.num_vox;
    %Collect data for brain region AIT    
    acc3(i) = x(i).apAnalyses(3).apROI{9}.acc;
    pval3(i) = x(i).apAnalyses(3).apROI{9}.pval;
    perm_upper3(i) = x(i).apAnalyses(3).apROI{9}.perm_upper;
    num_obs3(i) = x(i).apAnalyses(3).apROI{9}.num_obs;
    num_vox3(i) = x(i).apAnalyses(3).apROI{9}.num_vox;
end 

% %Plot Classifier Accuracy for each Brain Region
% plot(age,acc1+0.01,'-xb','linewidth',5); hold on;
% plot(age,acc2,'-xr','linewidth',5);
% plot(age,acc3-0.01,'-xg','linewidth',5);
% plot(age,perm_upper1,'--xb','linewidth',3);
% plot(age,perm_upper2,'--xr','linewidth',3);
% plot(age,perm_upper3,'--xg','linewidth',3);
% ylim([0,1.1]);
% xlabel('Age (Months)','fontsize',22,'interpreter','latex');
% ylabel('Accuracy Rate','fontsize',24,'interpreter','latex');
% tl = 'Classifier Accuracy by Brain Region';
% filename = 'Brain_Regions_Acc.png';
% title(tl,'fontsize',24,'interpreter','latex')
% h=legend('V1 Rate','PIT Rate','AIT Rate','V1 Threshold','PIT Threshold','AIT Threshold');set(h,'interpreter','latex','fontsize',12,'Location','southeast')
% set(gca,'fontsize',16);
% saveas(gcf, filename);
% close all;

all_accs = [acc1, acc2, acc3];
all_obs = [num_obs1, num_obs2, num_obs3];
all_vox = [num_vox1, num_vox2, num_vox3];

%Plot Classifier Accuracy by Number of Observations
plot(all_obs,all_accs,'xb','linewidth',5);
ylim([0.8,1.01]);
xlabel('Number of Observations','fontsize',22,'interpreter','latex');
ylabel('Accuracy Rate','fontsize',24,'interpreter','latex');
tl = 'Accuracy vs. Observations';
filename = 'Acc_vs_Obs.png';
title(tl,'fontsize',24,'interpreter','latex')
set(gca,'fontsize',16);
saveas(gcf, filename);
close all;

%Plot Classifier Accuracy by Number of Dimensions
plot(all_vox,all_accs,'xb','linewidth',5);
ylim([0.8,1.01]);
xlabel('Number of Dimensions','fontsize',22,'interpreter','latex');
ylabel('Accuracy Rate','fontsize',24,'interpreter','latex');
tl = 'Accuracy vs. Dimensions';
filename = 'Acc_vs_Dims.png';
title(tl,'fontsize',24,'interpreter','latex')
set(gca,'fontsize',16);
saveas(gcf, filename);
close all;