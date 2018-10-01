% Comparison plots for ananalytical vs. ABC banana

addpath('../dram')
addpath('../dram/utils')
addpath('../abc')

%% M-H banana
plot_test = bananatestFunction(1);

[bandwidth,density,X,Y] = kde2d(plot_test,2^8,[-5,-5],[5,15]);

%figure(1)
% subplot(1,2,1)
% scatter(plot_test(1:1000,1),plot_test(1:1000,2),'k.')
% hold on 
% bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
% bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
% box on
% xlim([-5,5])
% ylim([-5,15])
% xlabel('X')
% ylabel('Y')
% text(4,-4,'a','FontSize',20)

figure(2)
%subplot(1,2,2)
colormap(brewermap([],'Blues'))
imagesc(X(1,:)',Y(:,1),density)
set(gca,'Ydir','normal')
xlabel('X')
ylabel('Y')
matlab2tikz('my_fig.tex','standalone',true)
%text(4,-4,'b','FontSize',20)

%set(gcf,'PaperType','uslegal','PaperOrientation','landscape')

%print -dpdf -painters -fillpage MH-banana.pdf

%% Analytical vs. MCMC comparison

% labels = {'a','b','c','d','e','f','g','h'};
% methods = {'MH','AM','DR','DRAM','ABC-MH','ABC-AM','ABC-DR','ABC-DRAM'};
% 
% figure(1)
% for ii = 1:4
%    
%     analytical = bananatestFunction(ii);
%     [bandwidth,density,X,Y] = kde2d(analytical,2^8,[-5,-5],[5,15]);
%     
%     subplot(4,2,ii)
%     
%     colormap(brewermap([],'Blues'))
%     imagesc(X(1,:)',Y(:,1),density)
%     hold on
%     set(gca,'Ydir','normal')
%     %xlabel('X')
%     %ylabel('Y')
%     bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
%     bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
%     text(4,-3,labels{ii},'FontSize',12)
%     text(-4,12.5,methods{ii},'FontSize',12)
%     xlim([-5,5])
%     ylim([-5,15])
%     
% end
% 
% for ii = 1:4
%    
%    abc = ABCBananaTestFunction(ii);
%    
%    subplot(4,2,4+ii)
%    
%    bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
%    hold on
%    bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
%    bananaplot(median(abc(:,1)),median(abc(:,2)),median(abc(:,3)),...
%        median(abc(:,4)),median(abc(:,5)),1,1,'k-',0.5)
%    bananaplot(median(abc(:,1)),median(abc(:,2)),median(abc(:,3)),...
%        median(abc(:,4)),median(abc(:,5)),1,1,'k-',0.95)
%    xlim([-5,5])
%    ylim([-5,15])
%     %xlabel('X')
%     %ylabel('Y')
%    box on
%    text(4,-3,labels{4+ii},'FontSize',12)
%    text(-4,12.5,methods{4+ii},'FontSize',12)
%     
%     
% end
% 
% print -dpdf -painters -fillpage Analytical-vs-ABC.pdf
