function [] = plotbivariate(X,Y,xlab,ylab)
% Function to for the plotting of two random variables against each other.
% The plot produced is a bivariate scatter graph with two marginal
% histograms for both X and Y.
%
% INPUT:
% X = random variable 1
% Y = random variable 2

figure
%% Joint distribution
ax1 = subplot(4,4,[2,3,4,6,7,8,10,11,12]);
scatter(X,Y,'k','marker','.')
box on
grid on

%% Marginal distribution Y
ax2 = subplot(4,4,[1,5,9]);
h1 = histogram(Y,'normalization','pdf','facecolor','k','facealpha',0.3);
xlabel(ylab)
hold on

set(gca,'view',[-90,90],'xlim',get(ax1,'ylim'),'xtick',[],'ytick',[],'fontsize',12)

%% Marginal distribution X
ax3 = subplot(4,4,[14,15,16]);
h2 = histogram(X,'normalization','pdf','facecolor','k','facealpha',0.3);
xlabel(xlab)
hold on

set(gca,'xtick',[],'ytick',[],'fontsize',12)
linkaxes([ax1,ax3],'x')