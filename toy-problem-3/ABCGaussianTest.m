% Test of using ABC-MCMC with the forward solver abclikelihoodGauss.m
% This script compares the performance of ABC-MCMC with a uniform kernel
% to ABC-MCMC with a Gaussian kernel
% The parameters to a Gaussian distribution are estimated
% Author: Tom Connell
% Date: March 2018

clear;clc

addpath('../abc/abcutils')
addpath('../dram')
addpath('../dram/utils')


%% Loop over inference for kernel type

% 1 = Gaussian, 2 = Uniform
for kerneltype = 1:2
    
    % Use regular MCMC sampler, hence
    % No adaptive MCMC
    adaptint = 0;
    % No delayed rejection
    drscale = 0;

    % Number of iterations
    nsimu = 1000000;

    % Chain starting location
    % mu-x, mu-y, var-x, var-y, correlation
    start = [5,5,5,5,-0.5];

    % Parameter Bounds [xmin,ymin;xmax,ymax]
    bounds = [0,0,0,0,-1;10,10,10,10,1];
    
    % Set ABC tolerance
    tolerance = [1;1;1;1;1];
    
    % The forward problem
    if kerneltype == 1 % Gaussian
        forward = @(par,data) abclikelihoodgauss(par,data,tolerance,'g');
    end
    
    if kerneltype == 2 % Uniform
        forward = @(par,data) abclikelihoodgauss(par,data,tolerance,'u');
    end
    
    % Convariance Matrix of Markov chain transition kernel (Gaussian)
    qcov = eye(5);
    qcov(5,5) = 0.2;

    % Define passed structs
    clear model data params options

    model.ssfun = forward;

    % define the data
    % target mu-x, mu-y, var-x, var-y, correlation
    data = [2.5;7.5;4;6;2.4495];

    params.par0 = start;
    params.bounds = bounds;

    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;

    % Call to dramrun
    if kerneltype == 1
        [gresults,gchain] = dramrun(model,data,params,options);
    end
    
    if kerneltype == 2
        [uresults,uchain] = dramrun(model,data,params,options);
    end
    
end


%% Plotting

% Loop to produce two similar plots, one of the first 10,000 steps,
% the other of all 1,000,000
for ii = 1:2
    
    % How far into the chain to plot?
    if ii == 1
        n_to_plot = 10000;
        interval = 10;
    end
    
    if ii == 2
        n_to_plot = 1000000;
        interval = 1000;
    end
    
    figure(1)

    colors = brewermap(4,'Set1');
    
    % subplot 1
    subplot(6,2,1)
    % from abc folder
    PlotKernelComparison()
    xlabel('distance')
    ylabel('support')
    % set(gca,'YTick',[],'XAxisLocation','Top')
    set(gca,'YTick',[])
    text(0.925,0.15,'a','Units','normalized')

    % subplot 2
    subplot(6,2,2)
    
    [AX95,AY95] = error_ellipse([2.5,7.5],[4,2.4495;2.4495,6],chiqf_m(0.95,2));
    [AX50,AY50] = error_ellipse([2.5,7.5],[4,2.4495;2.4495,6],chiqf_m(0.95,2));
    plot(AX95,AY95,'r--','linewidth',1)
    hold on
    plot(AX50,AY50,'r--','linewidth',1)

    gCovar = sqrt(mean(gchain(:,3))*mean(gchain(:,4)))*mean(gchain(:,5));
    gCM = [mean(gchain(:,3)),gCovar;gCovar,mean(gchain(:,4))];
    [GX95,GY95] = error_ellipse([mean(gchain(:,1)),mean(gchain(:,2))],...
        gCM,chiqf_m(0.95,2));
    [GX50,GY50] = error_ellipse([mean(gchain(:,1)),mean(gchain(:,2))],...
        gCM,chiqf_m(0.5,2));
    plot(GX95,GY95,'color',colors(2,:),'linewidth',1)
    plot(GX50,GY50,'color',colors(2,:),'linewidth',1)

    uCovar = sqrt(mean(uchain(:,3))*mean(uchain(:,4)))*mean(uchain(:,5));
    uCM = [mean(uchain(:,3)),uCovar;uCovar,mean(uchain(:,4))];
    [UX95,UY95] = error_ellipse([mean(uchain(:,1)),mean(uchain(:,2))],...
        gCM,chiqf_m(0.95,2));
    [UX50,UY50] = error_ellipse([mean(uchain(:,1)),mean(uchain(:,2))],...
        gCM,chiqf_m(0.5,2));
    plot(UX95,UY95,'color',colors(4,:),'linewidth',1)
    plot(UX50,UY50,'color',colors(4,:),'linewidth',1)
    xlabel('X')
    ylabel('Y')
    %set(gca,'XAxisLocation','Top')
    text(0.925,0.15,'b','Units','normalized')
    
    % subplot 3
    subplot(6,2,3)
    histogram(gchain(1:n_to_plot,1),20,'facecolor',colors(2,:),'facealpha',0.5,'normalization','pdf')
    hold on
    histogram(uchain(1:n_to_plot,1),20,'facecolor',colors(4,:),'facealpha',0.5,'normalization','pdf')
    set(gca,'YTick',[])
    xlabel('\mu_X')
    xlim([0,10])
    text(0.925,0.15,'c','Units','normalized')
    % text(0.05,0.85,'\mu_X','Units','normalized')

    % subplot 4
    subplot(6,2,4)
    plot(1:interval:n_to_plot,gchain(1:interval:n_to_plot,1),'color',colors(2,:))
    hold on
    plot(1:interval:n_to_plot,uchain(1:interval:n_to_plot,1),'color',colors(4,:))
    rightPar = ones(1,n_to_plot/interval)*2.5;
    plot(1:interval:n_to_plot,rightPar,'r--')
    ylim([0,10])
    ylabel('\mu_X')
    xlabel('time')
    set(gca,'XTick',[])
    text(0.925,0.15,'d','Units','normalized')

    % subplot 5
    subplot(6,2,5)
    histogram(gchain(1:n_to_plot,2),20,'facecolor',colors(2,:),'facealpha',0.5,'normalization','pdf')
    hold on
    histogram(uchain(1:n_to_plot,2),20,'facecolor',colors(4,:),'facealpha',0.5,'normalization','pdf')
    set(gca,'YTick',[])
    xlabel('\mu_Y')
    xlim([0,10])
    text(0.925,0.15,'e','Units','normalized')
    % text(0.05,0.85,'\mu_Y','Units','normalized')
    
    % subplot 6
    subplot(6,2,6)
    plot(1:interval:n_to_plot,gchain(1:interval:n_to_plot,2),'color',colors(2,:))
    hold on
    plot(1:interval:n_to_plot,uchain(1:interval:n_to_plot,2),'color',colors(4,:))
    rightPar = ones(1,n_to_plot/interval)*7.5;
    plot(1:interval:n_to_plot,rightPar,'r--')
    ylim([0,10])
    ylabel('\mu_Y')
    xlabel('time')
    set(gca,'XTick',[])
    text(0.925,0.15,'f','Units','normalized')

    % subplot 7
    subplot(6,2,7)
    histogram(gchain(1:n_to_plot,3),20,'facecolor',colors(2,:),'facealpha',0.5,'normalization','pdf')
    hold on
    histogram(uchain(1:n_to_plot,3),20,'facecolor',colors(4,:),'facealpha',0.5,'normalization','pdf')
    set(gca,'YTick',[])
    xlabel('\sigma^2_X')
    xlim([0,10])
    text(0.925,0.15,'g','Units','normalized')
    % text(0.05,0.85,'\sigma^2_X','Units','normalized')

    % subplot 8
    subplot(6,2,8)
    plot(1:interval:n_to_plot,gchain(1:interval:n_to_plot,3),'color',colors(2,:))
    hold on
    plot(1:interval:n_to_plot,uchain(1:interval:n_to_plot,3),'color',colors(4,:))
    rightPar = ones(1,n_to_plot/interval)*4;
    plot(1:interval:n_to_plot,rightPar,'r--')
    ylim([0,10])
    ylabel('\sigma^2_X')
    xlabel('time')
    set(gca,'XTick',[])
    text(0.925,0.15,'h','Units','normalized')

    % subplot 9
    subplot(6,2,9)
    histogram(gchain(1:n_to_plot,4),20,'facecolor',colors(2,:),'facealpha',0.5,'normalization','pdf')
    hold on
    histogram(uchain(1:n_to_plot,4),20,'facecolor',colors(4,:),'facealpha',0.5,'normalization','pdf')
    set(gca,'YTick',[])
    xlabel('\sigma^2_Y')
    xlim([0,10])
    text(0.925,0.15,'i','Units','normalized')
    % text(0.05,0.85,'\sigma^2_Y','Units','normalized')

    % subplot 10
    subplot(6,2,10)
    plot(1:interval:n_to_plot,gchain(1:interval:n_to_plot,4),'color',colors(2,:))
    hold on
    plot(1:interval:n_to_plot,uchain(1:interval:n_to_plot,4),'color',colors(4,:))
    rightPar = ones(1,n_to_plot/interval)*6;
    plot(1:interval:n_to_plot,rightPar,'r--')
    ylim([0,10])
    ylabel('\sigma^2_Y')
    xlabel('time')
    set(gca,'XTick',[])
    text(0.925,0.15,'j','Units','normalized')

    % subplot 11
    subplot(6,2,11)
    histogram(gchain(1:n_to_plot,5),20,'facecolor',colors(2,:),'facealpha',0.5,'normalization','pdf')
    hold on
    histogram(uchain(1:n_to_plot,5),20,'facecolor',colors(4,:),'facealpha',0.5,'normalization','pdf')
    set(gca,'YTick',[])
    xlabel('\rho')
    text(0.925,0.15,'k','Units','normalized')

    % subplot 12
    subplot(6,2,12)
    plot(1:interval:n_to_plot,gchain(1:interval:n_to_plot,5),'color',colors(2,:))
    hold on
    plot(1:interval:n_to_plot,uchain(1:interval:n_to_plot,5),'color',colors(4,:))
    rightPar = ones(1,n_to_plot/interval)*0.5;
    plot(1:interval:n_to_plot,rightPar,'r--')
    ylim([-1,1])
    ylabel('\rho')
    xlabel('time')
    set(gca,'XTick',[])
    text(0.925,0.15,'l','Units','normalized')
    
    if ii == 1
        
        set(gcf,'units','centimeters','position',[0,0,15,20],'papersize',[15,20])
        print -dpdf -painters initalchains.pdf
        
    end
    
    if ii == 2
        
        set(gcf,'units','centimeters','position',[0,0,15,20],'papersize',[15,20])
        print -dpdf -painters fullschains.pdf    
    
    end
    
    % matlab2tikz(['gauss',num2str(ii),'.tex'],'width','\fwidth','height','\fheight')
    clf

end

% Close up demonstration of the initialization problem
figure(2)
CM = [4,2.4495;2.4495,6];
R = chol(CM);
obs = randn(100,2)*R;

subplot(1,2,1)
scatter(obs(:,1)+2.5,obs(:,2)+7.5,'k.')
hold on
[X95,Y95] = error_ellipse([2.5;7.5],CM,chiqf_m(0.95,2));
plot(X95,Y95,'r--')
[X50,Y50] = error_ellipse([2.5,7.5],CM,chiqf_m(0.5,2));
plot(X50,Y50,'r--')
box on
legend({'Observations','Causitive model'},'location','northwest')
legend boxoff
xlabel('X')
ylabel('Y')
text(0.9,0.1,'a','Units','normalized')

subplot(1,2,2)
plot(1:10:10000,uchain(1:10:10000,1),'k')
hold on
plot(1:10:10000,uchain(1:10:10000,2),'color',[0.7,0,0])
plot(1:10:10000,uchain(1:10:10000,3),'color',[0.0,0.5,0])
plot(1:10:10000,uchain(1:10:10000,4),'color',[0,0,0.7])
hold off
xlim([1,10000])
ylim([0,10])
legend({'\mu_X','\mu_Y','\sigma_X','\sigma_Y'},'location','northwest')
legend boxoff
xlabel('time')
ylabel('parameter value')
text(0.9,0.1,'b','Units','normalized')
%set(gca,'XTick',[])

set(gcf,'units','centimeters','position',[0,0,18.75,7.5],'papersize',[18.75,7.5])
print -dpdf -painters initqualms.pdf   

% matlab2tikz('init_qualms.tex','width','\fwidth','height','\fheight')