function [] = plotGaussVsUnif(chain)

figure

%% Plot data and bananafit
data = csvread('bananaData.dat');

ax1 = subplot(4,7,[4,5,6,7,11,12,13,14,18,19,20,21]);
scatter(data(:,1),data(:,2),'markerfacecolor',[0,0,1],'marker','.')
hold on
bananaplot(mean(chain(:,1)),mean(chain(:,2)),mean(chain(:,3)),...
    mean(chain(:,4)),mean(chain(:,5)),1,1,'k-')
xlabel('X')
ylabel('Y')
set(gca,'xtick',[],'ytick',[],'XAxisLocation', 'top','YAxisLocation', 'right')
legend({'Observed data','Posterior mean'})
legend boxoff
box on
grid on

%% muX

ax2 = subplot(4,7,1);
h1 = histogram(chain(:,1),'normalization','pdf','facecolor','k','facealpha',0.3);
hold on
set(gca,'view',[-90,90],'xlim',[-5,5],'xtick',[],'ytick',[],'fontsize',10)
xlabel('\mu_x')

ax22 = subplot(4,7,[2,3]);
plot(1:1:length(chain),chain(:,1),'k-')
set(gca,'ylim',[-5,5])

%% muY
ax3 = subplot(4,7,8);
h1 = histogram(chain(:,2),'normalization','pdf','facecolor','k','facealpha',0.3);
hold on
set(gca,'view',[-90,90],'xlim',[-5,5],'xtick',[],'ytick',[],'fontsize',10)
xlabel('\mu_Y')

ax33 = subplot(4,7,[9,10]);
plot(1:1:length(chain),chain(:,2),'k-')
set(gca,'ylim',[-5,5])

%% sigmaX

ax4 = subplot(4,7,15);
h1 = histogram(chain(:,3),'normalization','pdf','facecolor','k','facealpha',0.3);
hold on
set(gca,'view',[-90,90],'xlim',[0,10],'xtick',[],'ytick',[],'fontsize',10)
xlabel('\sigma^2_X')

ax44 = subplot(4,7,[16,17]);
plot(1:1:length(chain),chain(:,3),'k-')
set(gca,'ylim',[0,10])

%% simgaY

ax5 = subplot(4,7,22);
h1 = histogram(chain(:,4),'normalization','pdf','facecolor','k','facealpha',0.3);
hold on
set(gca,'view',[-90,90],'xlim',[0,10],'xtick',[],'ytick',[],'fontsize',10)
xlabel('\sigma^2_Y')

ax55 = subplot(4,7,[23,24]);
plot(1:1:length(chain),chain(:,4),'k-')
set(gca,'ylim',[0,10])

%% Correlation

ax6 = subplot(4,7,25);
h1 = histogram(chain(:,5),'normalization','pdf','facecolor','k','facealpha',0.3);
hold on
set(gca,'view',[-90,90],'xlim',[0,1],'xtick',[],'ytick',[],'fontsize',10)
xlabel('\rho')

ax66 = subplot(4,7,[26,27]);
plot(1:1:length(chain),chain(:,5),'k-')
set(gca,'ylim',[0,1])










