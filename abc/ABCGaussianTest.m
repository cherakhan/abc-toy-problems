% Test of using DRAM with the forward solver abclikelihoodGauss.m. 
% Author: Tom Connell
% Date: March 2018

clear;clc

addpath('./abcutils')
addpath('../dram')
addpath('../dram/utils')

% for all methods
methods = {'mh','am','dr','dram'};

%% Run all methods
for iMethods = 1:4
    
    method = methods{iMethods};
    
    switch method
        case 'mh'
            adaptint = 0;
            drscale = 0;
            plotTitle = 'M-H';
        case 'am'
            adaptint = 100;
            drscale = 0;
            plotTitle = 'AM';
        case 'dr'
            adaptint = 0;
            drscale = 5;
            plotTitle = 'DR';
        case 'dram'
            adaptint = 100;
            drscale = 5;
            plotTitle = 'DRAM';
    end
            
    %% Define input variables

    % Number of iterations
    nsimu = 50000;

    % Chain starting location
    start = [5,5,5,5,-0.5];

    % bounds [xmin,ymin;xmax,ymax]
    bounds = [-Inf,-Inf,0,0,-1;Inf,Inf,Inf,Inf,1];
    
    tolerance = [1;1;1;1;1];
    % The forward problem - ENCODES TOLERANCE
    forward = @(par,data) abclikelihoodgauss(par,data,tolerance);

    qcov = eye(5);
    qcov(5,5) = 0.2;

    %% Define passed structs

    clear model data params options

    model.ssfun = forward;

    % define the data
    % target muX,muY,varX,varY,correlation
    data = [2.5;7.5;4;6;2.4495];

    params.par0 = start;
    params.bounds = bounds;

    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;


    %% Call to dramrun

    [results,chain] = dramrun(model,data,params,options);
    
    %% Output MCMC chain results
    
    %csvwrite('/home/tom/thesis/dram/abc/output/guassian5par/gaussian/dramchain2.csv',chain)
    
    %% Plotting
    figure(1);
    
    subplot(4,3,(3*iMethods)-2);
    %subplot(1,3,1);
    scatter(chain(:,1),chain(:,2),'k','marker','.')
    box on
    grid on
    title(plotTitle)
    xlabel('\mu_x')
    ylabel('\mu_y')
    xlim([0,10])
    ylim([0,10])
    text(0.5,1,sprintf('\\tau = %4.1f',mean(iact(chain))))
    
    subplot(4,3,(3*iMethods)-1)
    %subplot(1,3,2)
    scatter(chain(:,3),chain(:,4),'k','marker','.')
    box on
    grid on
    title(plotTitle)
    xlabel('\sigma^2_x')
    ylabel('\sigma^2_y')
    xlim([0,10])
    ylim([0,10])
    text(0.5,1,sprintf('Accepted %3.1f %%',(results.accepted)*100))
    
    subplot(4,3,(3*iMethods));
    %subplot(1,3,3)
    plot(1:1:nsimu,chain(:,1),'k')
    hold on
    plot(1:1:nsimu,chain(:,2),'color',[0.7,0,0])
    plot(1:1:nsimu,chain(:,3),'color',[0.0,0.5,0])
    plot(1:1:nsimu,chain(:,4),'color',[0,0,0.7])
    hold off
    title(plotTitle)
    xlim([1,nsimu])
    ylim([0,10])
end
