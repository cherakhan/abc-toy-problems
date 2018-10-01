function [passed_chain] = ABCBananaTestFunction(method_to_use)

% Test DRAM with banana ABC problem
% Author: Tom Connell
% Date: April 2018

addpath('../dram')
addpath('../dram/utils')
addpath('../abc/abcutils')

% For all methods

methods = {'ABC-MH','ABC-AM','ABC-DR','ABC-DRAM'};


%% Run all the methods?

for iMethods = method_to_use:method_to_use
    
    %% Number of interations
    
    nsimu = 900000;
    
    
    %% Sampling method
    
    method = methods{iMethods};
    
    switch method
        case 'ABC-MH'
            adaptint = 0;
            drscale = 0;
        case 'ABC-AM'
            adaptint = nsimu/6;
            drscale = 0;
        case 'ABC-DR'
            adaptint = 0;
            drscale = 3;
        case 'ABC-DRAM'
            adaptint = nsimu/6;
            drscale = 3;
    end
    
    
    %% Means only - 2 par problem
    
%     % Starting location
%     start = [-2,2];
%     
%     % Parameter range
%     bounds = [-5,-5;5,5];
%     
%     % The likelihood approximation
%     tolerance = [0.1,0.1];    
%     bananaabc = @(parameters,data) abcbananafun([parameters,1,1,0.9,1,1],data,tolerance);
%     
%     % Transition kernel
%     qcov = eye(2)*0.1;
    
    
    %% Gaussian parameters only - 5 par problem
    
    % Starting location
    start = [1,2,2,2,0.6];
    
    % Parameter range
    bounds = [-5,-5,0.1,0.1,0;5,5,10,10,1];
    
    % The likelihood approximation
    normTerm = csvread('normTermGaussPars.csv');
    tolerance = [0.005,0.025,0.025,0.025,0.025,0.025,0.025,0.025];
    tolerance = tolerance*20;
    %tolerance = tolerance*5;
    tolerance = tolerance.*normTerm;
    bananaabc = @(parameters,data) abcbananafun([parameters,1,1],data,tolerance);
    
    % Transition kernel
    qcov = eye(5)*2;
    qcov(5,5) = 0.2;
%     qcov = csvread('adaptedqcov.csv');
    
    %% All paramters - 7 par problem
    
%     % Starting location
%     start = [1,2,2,2,0,0.5,1];
%     
%     % Parameter range
%     bounds = [-5,-5,0,0,-1,0,0;5,5,10,10,1,5,5];
%     
%     % The likelihood approximation
%     normTerm = csvread('normTermAllPars.csv');
%     tolerance = [0.005/5,0.025/5,0.025/5,0.025/5,0.025/5,0.025/5,0.025/5,0.025/5];
%     tolerance = tolerance.*normTerm;
%     bananaabc = @(parameters,data) abcbananafun(parameters,data,tolerance);
%     
%     % Transition kernel
%     qcov = eye(7)*0.1;
%     qcov(5,5) = 0.01;
    
    

    %% Define passed structs

    clear model data params options

    model.ssfun = bananaabc;

    % data to be passed to the forward
    data = csvread('bananaData.dat');

    params.par0 = start;
    params.bounds = bounds;

    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;
    
    
    %% Call to dramrun

    [results,chain] = dramrun(model,data,params,options);
    %[results,chain] = blockeddramrun(model,data,params,options);
    
    
    %% Plotting 
    % a number of routines for plotting different outputs
    
    
    %% 5 Parameter marginal posteriors + chain history
%     
%     figure(1)
%     
%     subplot(5,2,1)
%     histogram(chain(:,1))
%     subplot(5,2,2)
%     plot(1:1:nsimu,chain(:,1),'k')
%     
%     subplot(5,2,3)
%     histogram(chain(:,2))
%     subplot(5,2,4)
%     plot(1:1:nsimu,chain(:,2),'k')
%     
%     subplot(5,2,5)
%     histogram(chain(:,3))
%     subplot(5,2,6)
%     plot(1:1:nsimu,chain(:,3),'k')
%     
%     subplot(5,2,7)
%     histogram(chain(:,4))
%     subplot(5,2,8)
%     plot(1:1:nsimu,chain(:,4),'k')
%     
%     subplot(5,2,9)
%     histogram(chain(:,5))
%     subplot(5,2,10)
%     plot(1:1:nsimu,chain(:,5),'k')
    
    %% 5 parameter marginals + banana for thesis
    
%     figure(1)
%     
%     subplot(3,2,1)
%     histogram(chain(:,1),'normalization','pdf','facecolor','k','facealpha',0.5)
%     hold on
%     plot(linspace(0,0),linspace(0,3),'r-')
%     set(gca,'ytick',[])
%     xlabel('\mu_X')
%     
%     subplot(3,2,2)
%     histogram(chain(:,2),'normalization','pdf','facecolor','k','facealpha',0.5)
%     hold on
%     plot(linspace(0,0),linspace(0,0.8),'r-')
%     set(gca,'ytick',[])
%     xlabel('\mu_Y')
%     
%     subplot(3,2,3)
%     histogram(chain(:,3),'normalization','pdf','facecolor','k','facealpha',0.5)
%     hold on
%     plot(linspace(1,1),linspace(0,3),'r-')
%     set(gca,'ytick',[])
%     ylim([0,3])
%     xlabel('\sigma_X')
%     
%     subplot(3,2,4)
%     histogram(chain(:,4),'normalization','pdf','facecolor','k','facealpha',0.5)
%     hold on
%     plot(linspace(1,1),linspace(0,1),'r-')
%     set(gca,'ytick',[])
%     xlabel('\sigma_Y')
%     
%     subplot(3,2,5)
%     histogram(chain(:,5),'normalization','pdf','facecolor','k','facealpha',0.5)
%     hold on
%     plot(linspace(0.9,0.9),linspace(0,8),'r-')
%     set(gca,'ytick',[])
%     ylim([0,8])
%     xlabel('\rho')
%     
%     subplot(3,2,6)
%     bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
%     hold on
%     bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
%     bananaplot(median(chain(:,1)),median(chain(:,2)),...
%         median(chain(:,3)),median(chain(:,4)),...
%         median(chain(:,5)),1,1,'k-',0.5)
%     bananaplot(median(chain(:,1)),median(chain(:,2)),...
%         median(chain(:,3)),median(chain(:,4)),...
%         median(chain(:,5)),1,1,'k-',0.95)
%     xlim([-5,5])
%     ylim([-5,15])
%     xlabel('X')
%     ylabel('Y')
%     box on
    
    
    %% 7 Parameter maginal posteriors + chain history
    
%     figure(1)
%     
%     subplot(7,2,1)
%     histogram(chain(:,1))
%     subplot(7,2,2)
%     plot(1:1:nsimu,chain(:,1),'k')
%     
%     subplot(7,2,3)
%     histogram(chain(:,2))
%     subplot(7,2,4)
%     plot(1:1:nsimu,chain(:,2),'k')
%     
%     subplot(7,2,5)
%     histogram(chain(:,3))
%     subplot(7,2,6)
%     plot(1:1:nsimu,chain(:,3),'k')
%     
%     subplot(7,2,7)
%     histogram(chain(:,4))
%     subplot(7,2,8)
%     plot(1:1:nsimu,chain(:,4),'k')
%     
%     subplot(7,2,9)
%     histogram(chain(:,5))
%     subplot(7,2,10)
%     plot(1:1:nsimu,chain(:,5),'k')
%     
%     subplot(7,2,11)
%     histogram(chain(:,6))
%     subplot(7,2,12)
%     plot(1:1:nsimu,chain(:,6),'k')
%     
%     subplot(7,2,13)
%     histogram(chain(:,7))
%     subplot(7,2,14)
%     plot(1:1:nsimu,chain(:,7),'k')


    %% 5 Parameter marginal posteriors + bananaplot
    
%     figure(2)
%     
%     subplot(3,2,1)
%     histogram(chain(:,1),'Normalization','pdf')
%     hold on
%     plot(linspace(0,0),linspace(0,3),'k--','linewidth',2)
%     xlabel('\mu_X')
%     set(gca,'ytick',[])
%     
%     subplot(3,2,2)
%     histogram(chain(:,2),'Normalization','pdf')
%     hold on
%     plot(linspace(0,0),linspace(0,0.6),'k--','linewidth',2)    
%     xlabel('\mu_Y')
%     set(gca,'ytick',[])
%     
%     subplot(3,2,3)
%     histogram(chain(:,3),'Normalization','pdf')
%     hold on
%     plot(linspace(1,1),linspace(0,3),'k--','linewidth',2) 
%     xlabel('\sigma^2_X')
%     set(gca,'ytick',[])
% 
%     subplot(3,2,4)
%     histogram(chain(:,4),'Normalization','pdf')
%     hold on
%     plot(linspace(1,1),linspace(0,1.5),'k--','linewidth',2) 
%     xlabel('\sigma^2_Y')
%     set(gca,'ytick',[])
%     
%     subplot(3,2,5)
%     histogram(chain(:,5),'Normalization','pdf')
%     hold on
%     plot(linspace(0.9,0.9),linspace(0,6),'k--','linewidth',2) 
%     xlabel('\rho')
%     set(gca,'ytick',[])
%     text(0,5,['Acceptance Rate = ',num2str(results.accepted*100),'%'])
%        
%     subplot(3,2,6)
%     scatter(data(:,1),data(:,2),'k.')
%     hold on
%     bananaplot(0,0,1,1,0.9,1,1,'k-')
%     bananaplot(median(chain(:,1)),median(chain(:,2)),median(chain(:,3)),...
%         median(chain(:,4)),median(chain(:,5)),1,1,'r-')
%     xlabel('X')
%     ylabel('Y')
%     box on
%     grid on


    %% 7 Parameter marginal posteriors + bananaplot
    
%     figure(2)
%     
%     subplot(4,2,1)
%     histogram(chain(:,1),'Normalization','pdf')
%     hold on
%     plot(linspace(0,0),linspace(0,3),'k--','linewidth',2)
%     xlabel('\mu_X')
%     set(gca,'ytick',[])
%     
%     subplot(4,2,2)
%     histogram(chain(:,2),'Normalization','pdf')
%     hold on
%     plot(linspace(0,0),linspace(0,0.6),'k--','linewidth',2)    
%     xlabel('\mu_Y')
%     set(gca,'ytick',[])
%     
%     subplot(4,2,3)
%     histogram(chain(:,3),'Normalization','pdf')
%     hold on
%     plot(linspace(1,1),linspace(0,3),'k--','linewidth',2) 
%     xlabel('\sigma^2_X')
%     set(gca,'ytick',[])
% 
%     subplot(4,2,4)
%     histogram(chain(:,4),'Normalization','pdf')
%     hold on
%     plot(linspace(1,1),linspace(0,1.5),'k--','linewidth',2) 
%     xlabel('\sigma^2_Y')
%     set(gca,'ytick',[])
%     
%     subplot(4,2,5)
%     histogram(chain(:,5),'Normalization','pdf')
%     hold on
%     plot(linspace(0.9,0.9),linspace(0,6),'k--','linewidth',2) 
%     xlabel('\rho')
%     set(gca,'ytick',[])
%     text(0,5,['Acceptance Rate = ',num2str(results.accepted*100),'%'])
%     
%     subplot(4,2,6)
%     histogram(chain(:,6),'Normalization','pdf')
%     hold on
%     plot(linspace(1,1),linspace(0,1),'k--','linewidth',2) 
%     xlabel('Banana 1')
%     set(gca,'ytick',[])
%     
%     subplot(4,2,7)
%     histogram(chain(:,7),'Normalization','pdf')
%     hold on
%     plot(linspace(1,1),linspace(0,1),'k--','linewidth',2) 
%     xlabel('Banana 2')
%     set(gca,'ytick',[])
%     
%     subplot(4,2,8)
%     scatter(data(:,1),data(:,2),'k.')
%     hold on
%     bananaplot(0,0,1,1,0.9,1,1,'k-')
%     bananaplot(mean(chain(:,1)),mean(chain(:,2)),mean(chain(:,3)),...
%         mean(chain(:,4)),mean(chain(:,5)),mean(chain(:,6)),mean(chain(:,6)),'r-')
%     xlabel('X')
%     ylabel('Y')
%     box on
%     grid on
    
%% MCMC comparison plot
% figure(3)
% subplot(4,2,4+iMethods)
% %scatter(data(:,1),data(:,2),'r.')
% hold on
% bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
% bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
% bananaplot(median(chain(:,1)),median(chain(:,2)),median(chain(:,3)),...
%     median(chain(:,4)),median(chain(:,5)),1,1,'k-',0.5)
% bananaplot(median(chain(:,1)),median(chain(:,2)),median(chain(:,3)),...
%     median(chain(:,4)),median(chain(:,5)),1,1,'k-',0.95)
% xlim([-5,5])
% ylim([-5,15])
% %xlabel('X')
% %ylabel('Y')
% box on
% %grid on
% text(-4.5,-3,methods{iMethods})

%% AM comparison plot

% if strcmp(method,'ABC-AM')
%     
%     c1std = chiqf_m(0.95,2);
%     
%     figure(2)
%     
%     subplot(1,2,1)
%     plot(chain(:,1),chain(:,2),'k.')
%     hold on
%     
%     mu1 = [mean(chain(:,1)),mean(chain(:,2))];
%     CMua1 = [2,0;0,2;];
%     CMa1 = results.qcov(1:2,1:2);
%     
%     [xua1,yua1] = ellipse(mu1,CMua1*c1std);
%     [xa1,ya1] = ellipse(mu1,CMa1*c1std);
%     
%     plot(xua1,yua1,'r-','linewidth',1.5)
%     %hold on
%     plot(xa1,ya1,'b-','linewidth',1.5)
%     xlabel('\mu_X')
%     ylabel('\mu_Y')
%     ylim([-5,5])
%     xlim([-5,5])
%     
%     text(4,-4,'a','FontSize',20)
%     
%     subplot(1,2,2)
%     plot(chain(:,3),chain(:,4),'k.')
%     hold on
%     
%     mu2 = [median(chain(:,3)),median(chain(:,4))];
%     CMua2 = [2,0;0,2;];
%     CMa2 = results.qcov(3:4,3:4);
%     
%     [xua2,yua2] = ellipse(mu2,CMua2*c1std);
%     [xa2,ya2] = ellipse(mu2,CMa2*c1std);
%     
%     plot(xua2,yua2,'r-','linewidth',1.5)
%     %hold on
%     plot(xa2,ya2,'b-','linewidth',1.5)
%     
%     xlabel('\sigma^2_X')
%     ylabel('\sigma^2_Y')
%     
%     xlim([0,10])
%     ylim([0,10])
%     
%     text(9,1,'b','FontSize',20)
%     
%     legend({'Posterior samples','Initial proposal','Adapted proposal'})
%     legend boxoff
%     
%     print -dpdf -painters -fillpage Am-demonstration.pdf
%     
%     
% end


end

passed_chain = chain;