% This is a script to run ABC-MCMC compared to MCMC w/ likelihood for a
% linear regression problem. Data is in the data.dat file. The true
% solution is m = 5, b = 10, and sigma = 5.

clear;clc

addpath('~/thesis/dram/dram')
addpath('~/thesis/dram/dram/utils')
addpath('~/thesis/dram/abc/abcutils')
addpath('~/thesis/dram/abc/LinearRegression')

% Define the possible methods
methods = {'mh','am','dr','dram'};

%% ABC-MCMC
for iMethod = 4:4
    
    method = methods{iMethod};
    
    switch method
        case 'mh'
            adaptint = 0;
            drscale = 0;
        case 'am'
            adaptint = 1000;
            drscale = 0;
        case 'dr'
            adaptint = 0;
            drscale = 3;
        case 'dram'
            adaptint = 1000;
            drscale = 3;
    end
    
    %% Define input vars
    
    % Number of iterations
    nsimu = 10000;
    
    % Chain starting location
    start = [10,5];
    
    % bounds [xmin,ymin;xmax,ymax]
    bounds = [0,-20;20,20];
    
    % ABC component
    normTerm = csvread('normilizationTerm.csv');
    tolerance = [0.015,0.015];
    tolerance = tolerance.*normTerm;
    
    forward = @(par,data) abclikelihood([par,5],tolerance,data);
    
    % Transition kernel
    qcov = eye(2)*0.05;

    %% define passed structs
    
    clear model data params options
    
    model.ssfun = forward;
    
    data = csvread('data.csv');
    
    params.par0 = start;
    params.bounds = bounds;

    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;
    
    
    %% Call to DRAM
    
    [abcResults,abcChain] = dramrun(model,data,params,options);
    
end


%% Analytical Likelihood MCMC
for iMethod = 4:4
    
    method = methods{iMethod};
    
    switch method
        case 'mh'
            adaptint = 0;
            drscale = 0;
        case 'am'
            adaptint = 1000;
            drscale = 0;
        case 'dr'
            adaptint = 0;
            drscale = 3;
        case 'dram'
            adaptint = 1000;
            drscale = 3;
    end
    
    %% Define input vars
    
    % Number of iterations
    nsimu = 10000;
    
    % Chain starting location
    start = [10,5];
    
    % bounds [xmin,ymin;xmax,ymax]
    bounds = [0,-20;20,20];
    
    forward = @(par,data) analyticallikelihood([par,5],data);
    
    % Transition kernel
    qcov = eye(2)*0.05;

    %% define passed structs
    
    clear model data params options
    
    model.ssfun = forward;
    
    data = csvread('data.csv');
    
    params.par0 = start;
    params.bounds = bounds;

    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;
    
    
    %% Call to DRAM
    
    [anaResults,anaChain] = dramrun(model,data,params,options);
    
end

%% Plotting

figure(1)

subplot(3,1,1)
histogram(anaChain(:,1),'facecolor',[0,0,0.7],'facealpha',0.5,'Normalization','pdf')
hold on
histogram(abcChain(:,1),'facecolor',[1,0,0],'facealpha',0.5,'Normalization','pdf')
plot(linspace(5,5),linspace(0,20),'k--','linewidth',2)
legend({'ABC-MCMC','Analytical-Likelihood','True parameter value'},'location','northeast')
legend boxoff
xlabel('slope')

subplot(3,1,2)
histogram(anaChain(:,2),'facecolor',[0,0,0.7],'facealpha',0.5,'Normalization','pdf')
hold on
histogram(abcChain(:,2),'facecolor',[1,0,0],'facealpha',0.5,'Normalization','pdf')
plot(linspace(10,10),linspace(0,1),'k--','linewidth',2)
legend({['Acceptance Rate = ',num2str(abcResults.accepted*100),'%'],...
    ['Acceptance Rate = ',num2str(anaResults.accepted*100),'%'],'True parameter value'},'location','northwest')
legend boxoff

xlabel('intercept')

subplot(3,1,3)
scatter((-50:1:50)',data,25,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
hold on
% abc solution
x = [(-50:1:50)',ones(101,1)];
betaABC = [mean(abcChain(:,1));mean(abcChain(:,2))];
y = x*betaABC;
plot(x(:,1),y,'color',[0.5,0,0],'linewidth',2)
box on
legend({'Observed Data','ABC solution'},'location','northwest')
legend boxoff
xlabel('X')
ylabel('Y')