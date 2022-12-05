% Rerun on existing .mat file, to make adjustments to plot.

%% Nominal plots: Think about making these subplots (left-right)!
% z1 vs time:
figure(1)
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1.5;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.5;
plotHarc(tNom,jNom,xNom(:,1),[],modificatorF,modificatorJ);
hold on
plot(deltaVecNom(3),deltaVecNom(1),'k.','MarkerSize', 14)
strDeltaNom = [num2str(deltaVecNom(3)),'s'];
hold on
text(deltaVecNom(3),deltaVecNom(1),strDeltaNom,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
hold off
axis([0 4 -10 60]);
xlabel('t','FontSize',16)
ylabel('z_1','FontSize',16)

% axes('Position',[0.3 0.8 0.2 0.1])
% box on
% hold on
% plotHarc(tNom,jNom,xNom(:,1),[],modificatorF,modificatorJ);
% hold on
% plot(deltaVecNom(3),deltaVecNom(1),'k.','MarkerSize', 14)
% strDeltaNom = [num2str(deltaVecNom(3)),'s'];
% hold on
% text(deltaVecNom(3),deltaVecNom(1),strDeltaNom,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
% hold off
% set(gca,'xtick',[0 2.5 5])
% set(gca,'ytick',[0 25 60])
% axis([0 5 -10 60])
% hold off
% 
% % Log scale L(z_1)-L^* vs time:
% subplot(1,2,2),semilogy(tNomA,lNomA,'LineWidth',1.5);
% axis([0 2000 10^(-36) 10^(4)]);
% set(gca,'ytick',[10^(-36) 10^(-32) 10^(-28) 10^(-24) 10^(-20) 10^(-16) 10^(-12) 10^(-8) 10^(-4) 10^(0) 10^4])
% legend({'c_{1,0}=6819.67593','c_{1,0}=6819.676','c_{1,0}=6819.68','c_{1,0}=6820'},'Location','northeast')
% xlabel('t','FontSize',16)
% ylabel('L(z_1)-L^*','FontSize',16)
% 
% axes('Position',[0.65 0.4 0.2 0.1])
% box on
% hold on
% semilogy(tNomA,lNomA,'LineWidth',1.5);
% hold off
% set(gca,'xtick',[0 1 2])
% set(gca,'ytick',[10^(-20) 10^(-10) 10^(0)])
% axis([0 2 10^(-20) 10^(5)])
% hold off

saveas(gcf,'Plots\PlotNominalLog','epsc')