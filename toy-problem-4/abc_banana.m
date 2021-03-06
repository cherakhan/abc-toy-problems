% Test ABC-DRAM with banana problem

clear;clc

addpath('../dram')
addpath('../dram/utils')
addpath('../abc/abcutils')


%% Define problem

problem_type = {'Gaussian','fullbanana'};

% Here we only consider solving for the 5 Gaussian parameters
problem = problem_type{1};


%% Define observed data based

% Number of observations
n_obs = 1000;

% True parameter values
true_mux = 0;
true_muy = 0;
true_corr = 0.9;
true_varx = 1;
true_vary = 1;
true_bananity = [1,1];

% Simulate observed data
true_Mu = [true_mux;true_muy];
true_covar = true_corr*sqrt(true_varx*true_vary);
true_Sigma = [true_varx,true_covar;true_covar,true_vary];
data = bananafun(randn(n_obs,2)*chol(true_Sigma),true_bananity,1);

% Plot observed data
figure(1)
plotbivariate(data(:,1),data(:,2),'X','Y')
set(gcf,'units','centimeters','position',[0,0,20,15],'papersize',[20,15])
print('-dpdf','-painters','observeddata.pdf')


%% Define common sampler variables for ABC

% Number of time steps in the chain
nsimu = 100002;

% For targeting Gaussian dist. parameter only 
if strcmp('Gaussian',problem)
   
    % Chain starting location
    start = [(rand*10)-5,(rand*10)-5,rand*10,rand*10,rand];
    
    % Parameter range
    bounds = [-5,-5,0.1,0.1,0;5,5,10,10,1];
    
    % The likelihood approximation
    range = normterm(problem,data);
    tolerance = ones(1,length(range))*0.625;
    % tighten the tolerance on the first polynomial term
    tolerance(1) = tolerance(1)*0.2;
    
    % Forward
    bananaabc = @(parameters,data) abcbananafun([parameters,1,1],data,tolerance);
    
    % Transition kernel
    qcov = eye(5);
    qcov(5,5) = 0.1;
     
end


% For also targeting the banana transformation parameters
if strcmp('fullbanana',problem)
   
    % Chain starting location
    start = [(rand*10)-5,(rand*10)-5,rand*10,rand*10,rand,rand*5,rand*5];
    
    % Parameter range
    bounds = [-5,-5,0,0,-1,0,0;5,5,10,10,1,5,5];
    
    % The likelihood approximation
    range = normterm(problem,data);
    tolerance = ones(1,length(range))*0.01;
    % tighten the tolerance on the first polynomial term
    tolerance(1) = tolerance(1)*0.2;
    
    % Forward
    bananaabc = @(parameters,data) abcbananafun([parameters,1,1],data,tolerance);
    
    % Transition kernel
    qcov = eye(7);
    qcov(5,5) = 0.1;   
    
end



%% Analytical sampling

% labels for large figure comparing ABC to analytical inference over 
% all methods
labels = {'a','b','c','d','e','f','g','h'};
analytical_methods = {'MH','AM','DR','DRAM'};

% Loop over analytical solution for all methods and establish comparison
% plot
for ii = 1:4
   
    chain = bananatestfunction(ii);
    
    [bandwidth,density,X,Y] = kde2d(chain,2^7,[-5,-5],[5,15]);
    
    % initial test plot with analytical MH
    if ii == 1
        
        figure
        subplot(1,2,1)
        scatter(chain(1:1000,1),chain(1:1000,2),'k.')
        hold on 
        bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
        bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
        box on
        xlim([-5,5])
        ylim([-5,15])
        xlabel('X')
        ylabel('Y')
        text(0.9,0.1,'a','units','normalized')
        
        subplot(1,2,2)
        colormap(brewermap([],'Blues'))
        imagesc(X(1,:)',Y(:,1),density)
        set(gca,'Ydir','normal')
        xlabel('X')
        ylabel('Y')
        text(0.9,0.1,'b','units','normalized')
        
        set(gcf,'units','centimeters','position',[0,0,25,10],'papersize',[25,10])
        print('-dpdf', '-painters', 'analyticaltestplot.pdf')
        
    end
    
    figure(10)
    subplot(4,2,ii)
    colormap(brewermap([],'Blues'))
    imagesc(X(1,:)',Y(:,1),density)
    hold on
    set(gca,'Ydir','normal')
    %xlabel('X')
    %ylabel('Y')
    bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
    bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
    text(4,-3,labels{ii})
    text(-4,12.5,analytical_methods{ii})
    xlim([-5,5])
    ylim([-5,15])
    
end

%% ABC for all methods

methods = {'ABC-MH','ABC-AM','ABC-DR','ABC-DRAM'};

acceptance_rates = zeros(1,5);

for i_methods = 1:4
    
    clear results chain
    
    % Sampling method
    method = methods{i_methods};
    
    switch method
        case 'ABC-MH'
            adaptint = 0;
            drscale = 0;
        case 'ABC-AM'
            adaptint = nsimu/6;
            drscale = 0;
        case 'ABC-DR'
            adaptint = 0;
            drscale = 2;
        case 'ABC-DRAM'
            adaptint = nsimu/6;
            drscale = 2;
    end
    
    % Define passed structs

    clear model params options

    model.ssfun = bananaabc;

    % data to be passed to the forward
    % data = csvread('bananaData.dat');

    params.par0 = start;
    params.bounds = bounds;

    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;
    
    % Call to blocked dramrun
    % [results,chain] = blockeddramrun(model,data,params,options);
    [results,chain] = dramrun(model,data,params,options);
    
    acceptance_rates(i_methods) = results.accepted;
    
    if strcmp(method,'ABC-MH')
        csvwrite('abc-mh.csv',chain);
    end
    
    % Marginals + banana (for thesis)
    if strcmp(problem,'Gaussian') && strcmp(method,'ABC-MH')
        
        figure    
        subplot(3,2,1)
        h1 = histogram(chain(:,1),'normalization','pdf','facecolor','k','facealpha',0.5);
        hold on
        topoff = max(h1.Values)+0.2*max(h1.Values);
        plot(linspace(0,0),linspace(0,topoff),'r-')
        set(gca,'ytick',[])
        xlabel('\mu_X')
        xlim([-5,5])
        ylim([0,topoff])
        text(0.9,0.9,'a','units','normalized')

        subplot(3,2,2)
        h2 = histogram(chain(:,2),'normalization','pdf','facecolor','k','facealpha',0.5);
        hold on
        topoff = max(h2.Values)+0.2*max(h2.Values);
        plot(linspace(0,0),linspace(0,topoff),'r-')
        set(gca,'ytick',[])
        xlabel('\mu_Y')
        xlim([-5,5])
        ylim([0,topoff])
        text(0.9,0.9,'b','units','normalized')

        subplot(3,2,3)
        h3 = histogram(chain(:,3),'normalization','pdf','facecolor','k','facealpha',0.5);
        hold on
        topoff = max(h3.Values)+0.2*max(h3.Values);
        plot(linspace(1,1),linspace(0,topoff),'r-')
        set(gca,'ytick',[])
        ylim([0,3])
        xlabel('\sigma_X')
        xlim([0,10])
        ylim([0,topoff])
        text(0.9,0.9,'c','units','normalized')

        subplot(3,2,4)
        h4 = histogram(chain(:,4),'normalization','pdf','facecolor','k','facealpha',0.5);
        hold on
        topoff = max(h4.Values)+0.2*max(h4.Values);
        plot(linspace(1,1),linspace(0,topoff),'r-')
        set(gca,'ytick',[])
        xlabel('\sigma_Y')
        xlim([0,10])
        ylim([0,topoff])
        text(0.9,0.9,'d','units','normalized')

        subplot(3,2,5)
        h5 = histogram(chain(:,5),'normalization','pdf','facecolor','k','facealpha',0.5);
        hold on
        topoff = max(h5.Values)+0.2*max(h5.Values);
        plot(linspace(0.9,0.9),linspace(0,topoff),'r-')
        set(gca,'ytick',[])
        ylim([0,8])
        xlabel('\rho')
        xlim([0,1])
        ylim([0,topoff])
        text(0.1,0.9,'e','units','normalized')

        subplot(3,2,6)
        bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
        hold on
        bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
        bananaplot(median(chain(:,1)),median(chain(:,2)),...
            median(chain(:,3)),median(chain(:,4)),...
            median(chain(:,5)),1,1,'k-',0.5)
        bananaplot(median(chain(:,1)),median(chain(:,2)),...
            median(chain(:,3)),median(chain(:,4)),...
            median(chain(:,5)),1,1,'k-',0.95)
        xlim([-5,5])
        ylim([-5,15])
        xlabel('X')
        ylabel('Y')
        box on
        text(0.9,0.9,'f','units','normalized')

        set(gcf,'units','centimeters','position',[0,0,15,20],'papersize',[15,20])
        print('-dpdf','-painters',['banana',method,'.pdf'])
        
    end
    
if strcmp(method,'ABC-AM')
    
    c1std = chiqf_m(0.95,2);
    
    figure
    
    subplot(1,2,1)
    plot(chain(1:5:length(chain),1),chain(1:5:length(chain),2),'k.')
    hold on
    
    mu1 = [mean(chain(:,1)),mean(chain(:,2))];
    CMua1 = [2,0;0,2;];
    CMa1 = results.qcov(1:2,1:2);
    
    [xua1,yua1] = ellipse(mu1,CMua1*c1std);
    [xa1,ya1] = ellipse(mu1,CMa1*c1std);
    
    plot(xua1,yua1,'r-','linewidth',1.5)
    plot(xa1,ya1,'b-','linewidth',1.5)
    xlabel('\mu_X')
    ylabel('\mu_Y')
    ylim([-5,5])
    xlim([-5,5])
    
    text(0.9,0.1,'a','units','normalized')
    
    normal_ar = acceptancerate(chain(1:(nsimu/6),:));
    
    subplot(1,2,2)
    plot(chain(1:5:length(chain),3),chain(1:5:length(chain),4),'k.')
    hold on
    
    mu2 = [median(chain(:,3)),median(chain(:,4))];
    CMua2 = [2,0;0,2;];
    CMa2 = results.qcov(3:4,3:4);
    
    [xua2,yua2] = ellipse(mu2,CMua2*c1std);
    [xa2,ya2] = ellipse(mu2,CMa2*c1std);
    
    plot(xua2,yua2,'r-','linewidth',1.5)
    plot(xa2,ya2,'b-','linewidth',1.5)
    
    xlabel('\sigma^2_X')
    ylabel('\sigma^2_Y')
    
    xlim([0,10])
    ylim([0,10])
    
    text(0.9,0.1,'b','units','normalized')
    
    am_ar = acceptancerate(chain((nsimu/6)*5:length(chain),:));
    
    text(0.4,0.6,'Adapted acceptance rate','units','normalized')
    text(0.5,0.55,['= ',num2str(round(am_ar*1000)/10),'%'],'units','normalized')
    
    text(0.4,0.7,'Initial acceptance rate','units','normalized')
    text(0.5,0.65,['= ',num2str(round(normal_ar*1000)/10),'%'],'units','normalized')
    
    legend({'Posterior samples','Initial proposal','Adapted proposal'})
    legend boxoff
    
    set(gcf,'units','centimeters','position',[0,0,25,10],'papersize',[25,10])
    
    print -dpdf -painters amdemonstration.pdf
    
end

% Large comparison plot continued...
figure(10)
subplot(4,2,4+i_methods)
bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
hold on
bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
bananaplot(median(chain(:,1)),median(chain(:,2)),median(chain(:,3)),...
   median(chain(:,4)),median(chain(:,5)),1,1,'k-',0.5)
bananaplot(median(chain(:,1)),median(chain(:,2)),median(chain(:,3)),...
   median(chain(:,4)),median(chain(:,5)),1,1,'k-',0.95)
xlim([-5,5])
ylim([-5,15])
%xlabel('X')
%ylabel('Y')
box on
text(4,-3,labels{4+i_methods})
text(-4,12.5,methods{i_methods})
if i_methods == 4
    set(gcf,'units','centimeters','position',[0,0,15,20],'papersize',[15,20])
    print('-dpdf','-painters','abcanacomparison.pdf')
end

end



%% Work out MH acceptance rate under single par update with same setting as
% large cycle through above
clear results chain
% Sampling method
method = 'ABC-MH';

switch method
    case 'ABC-MH'
        adaptint = 0;
        drscale = 0;
    case 'ABC-AM'
        adaptint = nsimu/6;
        drscale = 0;
    case 'ABC-DR'
        adaptint = 0;
        drscale = 2;
    case 'ABC-DRAM'
        adaptint = nsimu/6;
        drscale = 2;
end

% Define passed structs

clear model params options

model.ssfun = bananaabc;

params.par0 = start;
params.bounds = bounds;

options.nsimu = nsimu;
options.adaptint = adaptint;
options.drscale = drscale;
options.qcov = qcov;

[results,~] = blockeddramrun(model,data,params,options);
acceptance_rates(5) = results.accepted;    

    
 
%% Extra accurate dram, push inference to it's limits

clear results chain range tolerance bananaabc

% Reset the tolerance used
% The likelihood approximation
range = normterm(problem,data);
tolerance = ones(1,length(range))*0.075;
% tighten the tolerance on the first polynomial term
tolerance(1) = tolerance(1)*0.2;

% Forward
bananaabc = @(parameters,data) abcbananafun([parameters,1,1],data,tolerance);

% Sampling method
method = 'ABC-DRAM';

switch method
    case 'ABC-MH'
        adaptint = 0;
        drscale = 0;
    case 'ABC-AM'
        adaptint = nsimu/6;
        drscale = 0;
    case 'ABC-DR'
        adaptint = 0;
        drscale = 2;
    case 'ABC-DRAM'
        adaptint = nsimu/6;
        drscale = 2;
end

% Define passed structs

clear model params options

model.ssfun = bananaabc;

params.par0 = start;
params.bounds = bounds;

options.nsimu = nsimu;
options.adaptint = adaptint;
options.drscale = drscale;
options.qcov = qcov;

% Call to blocked dramrun
% [results,chain] = blockeddramrun(model,data,params,options);
[results,chain] = blockeddramrun(model,data,params,options);

figure    
subplot(3,2,1)
h1 = histogram(chain(:,1),'normalization','pdf','facecolor','k','facealpha',0.5);
hold on
topoff = max(h1.Values)+0.2*max(h1.Values);
plot(linspace(0,0),linspace(0,topoff),'r-')
set(gca,'ytick',[])
xlabel('\mu_X')
xlim([-5,5])
ylim([0,topoff])
text(0.9,0.9,'a','units','normalized')

subplot(3,2,2)
h2 = histogram(chain(:,2),'normalization','pdf','facecolor','k','facealpha',0.5);
hold on
topoff = max(h2.Values)+0.2*max(h2.Values);
plot(linspace(0,0),linspace(0,topoff),'r-')
set(gca,'ytick',[])
xlabel('\mu_Y')
xlim([-5,5])
ylim([0,topoff])
text(0.9,0.9,'b','units','normalized')

subplot(3,2,3)
h3 = histogram(chain(:,3),'normalization','pdf','facecolor','k','facealpha',0.5);
hold on
topoff = max(h3.Values)+0.2*max(h3.Values);
plot(linspace(1,1),linspace(0,topoff),'r-')
set(gca,'ytick',[])
ylim([0,3])
xlabel('\sigma_X')
xlim([0,10])
ylim([0,topoff])
text(0.9,0.9,'c','units','normalized')

subplot(3,2,4)
h4 = histogram(chain(:,4),'normalization','pdf','facecolor','k','facealpha',0.5);
hold on
topoff = max(h4.Values)+0.2*max(h4.Values);
plot(linspace(1,1),linspace(0,topoff),'r-')
set(gca,'ytick',[])
xlabel('\sigma_Y')
xlim([0,10])
ylim([0,topoff])
text(0.9,0.9,'d','units','normalized')

subplot(3,2,5)
h5 = histogram(chain(:,5),'normalization','pdf','facecolor','k','facealpha',0.5);
hold on
topoff = max(h5.Values)+0.2*max(h5.Values);
plot(linspace(0.9,0.9),linspace(0,topoff),'r-')
set(gca,'ytick',[])
ylim([0,8])
xlabel('\rho')
xlim([0,1])
ylim([0,topoff])
text(0.1,0.9,'e','units','normalized')

subplot(3,2,6)
bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
hold on
bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
bananaplot(median(chain(:,1)),median(chain(:,2)),...
    median(chain(:,3)),median(chain(:,4)),...
    median(chain(:,5)),1,1,'k-',0.5)
bananaplot(median(chain(:,1)),median(chain(:,2)),...
    median(chain(:,3)),median(chain(:,4)),...
    median(chain(:,5)),1,1,'k-',0.95)
xlim([-5,5])
ylim([-5,15])
xlabel('X')
ylabel('Y')
box on
text(0.9,0.9,'f','units','normalized')

set(gcf,'units','centimeters','position',[0,0,15,20],'papersize',[15,20])
print('-dpdf','-painters',['banana',method,'.pdf'])


%% Print table with data

info = {'ABC-MH';'ABC-AM';'ABC-DR';'ABC-DRAM';'Single parameter update'};
table(info,acceptance_rates')