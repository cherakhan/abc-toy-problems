% script to plot rejection results

tolerance2 = csvread('mhn10000t2');
tolerance5 = csvread('mhn10000t5');
tolerance7 = csvread('mhn10000t7');
tolerance10 = csvread('mhn10000t10');

colors = brewermap(4,'set1');
tolerances = {'\epsilon = 2.5','\epsilon = 5.0','\epsilon = 7.5','\epsilon = 10.0'};
ar = {'AR = 0.2%','AR = 1.3%','AR = 2.6%','AR = 3.9%'};
labels = {'\mu_X','\mu_Y','\sigma^2_X','\sigma^2_Y','\rho'};
analytical = [2.5;7.5;5;5;0.5];

range = [1,1,1,1,2];

figure
for ii = 1:5
    subplot(2,3,ii)
    
    if ii ~= 5
        [t2,xt2] = ksdensity(tolerance2(:,ii),'bandwidth',0.2);
        [t5,xt5] = ksdensity(tolerance5(:,ii),'bandwidth',0.2);
        [t7,xt7] = ksdensity(tolerance7(:,ii),'bandwidth',0.2);
        [t10,xt10] = ksdensity(tolerance10(:,ii),'bandwidth',0.2);
    else
        [t2,xt2] = ksdensity(tolerance2(:,ii),'bandwidth',0.04);
        [t5,xt5] = ksdensity(tolerance5(:,ii),'bandwidth',0.04);
        [t7,xt7] = ksdensity(tolerance7(:,ii),'bandwidth',0.04);
        [t10,xt10] = ksdensity(tolerance10(:,ii),'bandwidth',0.04);
    end
    
    plot(xt2,t2,'color',colors(1,:),'linewidth',1.5)
    hold on
    box on
    plot(xt5,t5,'color',colors(2,:),'linewidth',1.5)
    plot(xt7,t7,'color',colors(3,:),'linewidth',1.5)
    plot(xt10,t10,'color',colors(4,:),'linewidth',1.5)
    
    y = linspace(0,range(ii));
    x = ones(length(y),1);
    x = x*analytical(ii);
    plot(x,y,'k--','linewidth',2)
    
    if ii == 1
        legend(tolerances,'location','northeast')
        legend boxoff
    end
    
    if ii == 4
       legend(ar,'location','northeast')
       legend boxoff
    end
    
    xlabel(labels{ii})
    ylabel('Posterior density')
    
    set(gca,'ytick',[])
    
end

subplot(2,3,6)
% tolerance 2
mu2 = [mean(tolerance2(:,1)),mean(tolerance2(:,2))];
varx2 = mean(tolerance2(:,3));
vary2 = mean(tolerance2(:,4));
corr2 = mean(tolerance2(:,5));
covar2 = corr2*sqrt(varx2*vary2);
[Xt2,Yt2] = error_ellipse(mu2,[varx2,covar2;covar2,vary2],5.991);
plot(Xt2,Yt2,'color',colors(1,:),'linewidth',1.5)
hold on
% tolerance 5
mu5 = [mean(tolerance5(:,1)),mean(tolerance5(:,2))];
varx5 = mean(tolerance5(:,3));
vary5 = mean(tolerance5(:,4));
corr5 = mean(tolerance5(:,5));
covar5 = corr5*sqrt(varx5*vary5);
[Xt5,Yt5] = error_ellipse(mu5,[varx5,covar5;covar5,vary5],5.991);
plot(Xt5,Yt5,'color',colors(2,:),'linewidth',1.5)
% tolerance 7
mu7 = [mean(tolerance7(:,1)),mean(tolerance7(:,2))];
varx7 = mean(tolerance7(:,3));
vary7 = mean(tolerance7(:,4));
corr7 = mean(tolerance7(:,5));
covar7 = corr7*sqrt(varx7*vary7);
[Xt7,Yt7] = error_ellipse(mu7,[varx7,covar7;covar7,vary7],5.991);
plot(Xt7,Yt7,'color',colors(3,:),'linewidth',1.5)
% tolerane 10
mu10 = [mean(tolerance10(:,1)),mean(tolerance10(:,2))];
varx10 = mean(tolerance10(:,3));
vary10 = mean(tolerance10(:,4));
corr10 = mean(tolerance10(:,5));
covar10 = corr10*sqrt(varx10*vary10);
[Xt10,Yt10] = error_ellipse(mu10,[varx10,covar10;covar10,vary10],5.991);
plot(Xt10,Yt10,'color',colors(4,:),'linewidth',1.5)
% analytical
[X,Y] = error_ellipse([2.5,7.5],[5,2.5;2.5,5],5.991);
plot(X,Y,'k--','linewidth',2)
ylabel('Y')
xlabel('X')