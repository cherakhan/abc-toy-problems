% This is a script to run ABC-MCMC compared to MCMC w/ likelihood for a
% linear regression problem. The true solution is 
% m = 5, b = 10, and sigma = 5.
% In my thesis it is used to show that ABC-MCMC improves the efficiency of 
% statistical inference compared to a rejection scheme

clear;clc;clf

addpath('../dram/utils')
addpath('../dram')
addpath('../abc/abcutils')


%% Create the data

% Create noisy data from a linear model
% 101 samples from -50:50
x = [(-50:1:50)',ones(101,1)];

% True parameters [slope;intercept]
beta = [5;10];

% Clean data computation
y = x*beta;

% Noise 
noise = randn(length(y),1)*5;

% Add noise to data
data = y+noise;


%% ABC-MCMC    
% Define input vars

% Use regular MCMC sampler, hence
% No adaptive MCMC
adaptint = 0;
% No delayed rejection
drscale = 0;

% Number of iterations
nsimu = 20000;

% Chain starting location
start = [7,11];

% Parameter bounds [xmin,ymin;xmax,ymax]
bounds = [0,-20;20,20];

% ABC tolerance
tolerance = [2;2];

% The approximate likelihood term (p(S|S*,theta*))
abc_likelihood = @(par,data) abclikelihood([par,5],tolerance,data);

% Convariance Matrix of Markov chain transition kernel (Gaussian)
qcov = eye(2);

% Define the structs passed to the sampling algorithm
clear model params options

model.ssfun = abc_likelihood;

params.par0 = start;
params.bounds = bounds;

options.nsimu = nsimu;
options.adaptint = adaptint;
options.drscale = drscale;
options.qcov = qcov;

% Call to the MCMC sampling algorithm
[abc_results,abc_chain] = dramrun(model,data,params,options);
    

%% Analytical Likelihood MCMC
% Define input vars

% Use regular MCMC sampler, hence
% No adaptive MCMC
adaptint = 0;
% No delayed rejection
drscale = 0;

% Number of iterations
nsimu = 20000;

% Chain starting location
start = [7,11];

% Parameter bounds [xmin,ymin;xmax,ymax]
bounds = [0,-20;20,20];

% The analytical likelihood term (p(y|theta))
likelihood = @(par,data) analyticallikelihood([par,5],data);

% Convariance Matrix of Markov chain transition kernel (Gaussian)
qcov = eye(2);

% define passed structs
clear model params options

model.ssfun = likelihood;

params.par0 = start;
params.bounds = bounds;

options.nsimu = nsimu;
options.adaptint = adaptint;
options.drscale = drscale;
options.qcov = qcov;

% Call to the MCMC sampling algorithm
[anaResults,anaChain] = dramrun(model,data,params,options);


%% Plotting

figure(1)

subplot(3,2,1)
histogram(anaChain(:,1),300,'facecolor',[0,0,0.7],'facealpha',0.5,'Normalization','pdf')
hold on
histogram(abc_chain(:,1),100,'facecolor',[1,0,0],'facealpha',0.5,'Normalization','pdf')
plot(linspace(5,5),linspace(0,30),'k--','linewidth',2)
legend({'Analytical-Likelihood','ABC-MCMC','True parameter value'},'location','northeast')
legend boxoff
xlabel('slope')
ylabel('posterior density')
set(gca,'YTick',[])
%text(6.8,2.5,'a','FontSize',20)
text(0.9,0.1,'a','Units','normalized')

subplot(3,2,2)
plot(1:2:length(abc_chain),abc_chain(1:2:length(abc_chain),1),'k-')
xlabel('time')
ylabel('slope')
set(gca,'XTick',[])
text(0.9,0.1,'b','Units','normalized')

subplot(3,2,3)
histogram(anaChain(:,2),'facecolor',[0,0,0.7],'facealpha',0.5,'Normalization','pdf')
hold on
histogram(abc_chain(:,2),'facecolor',[1,0,0],'facealpha',0.5,'Normalization','pdf')
plot(linspace(10,10),linspace(0,2.5),'k--','linewidth',2)
legend({['Acceptance Rate = ',num2str(abc_results.accepted*100),'%'],...
    ['Acceptance Rate = ',num2str(anaResults.accepted*100),'%'],'True parameter value'},'location','northeast')
legend boxoff
xlabel('intercept')
ylabel('posterior density')
set(gca,'YTick',[])
text(0.9,0.1,'c','Units','normalized')

subplot(3,2,4)
plot(1:2:length(abc_chain),abc_chain(1:2:length(abc_chain),2),'k-')
xlabel('time')
ylabel('intercept')
set(gca,'XTick',[])
text(0.9,0.1,'d','Units','normalized')

subplot(3,2,[5,6])
scatter((-50:1:50)',data,25,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
hold on
% abc solution
x = [(-50:1:50)',ones(101,1)];
betaABC = [mean(abc_chain(:,1));mean(abc_chain(:,2))];
y = x*betaABC;
plot(x(:,1),y,'color',[0.7,0,0],'linewidth',2)
box on
legend({'Observed Data','Mean ABC posterior'},'location','northwest')
legend boxoff
xlabel('X')
ylabel('Y')
text(0.95,0.1,'e','Units','normalized')

set(gcf,'units','centimeters','position',[0,0,18.75,25],'papersize',[18.75,25])
print -dpdf -painters linearregression.pdf