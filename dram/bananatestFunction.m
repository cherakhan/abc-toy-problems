function [passed_chain] = bananatestFunction(method_to_use)

% Test DRAM/MCMC with a banana shaped distribution

addpath('./utils');
addpath('../abc')

% run all the methods
methods = {'MH','AM', 'DR', 'DRAM'};
for mi=method_to_use:method_to_use
  
  method = methods{mi};

  nsimu = 90000;
                                                          
  switch method
   case 'MH'
    drscale  = 0;
    adaptint = 0;
   case 'DR'
    drscale  = 2;
    adaptint = 0;
   case 'AM'
    drscale  = 0;
    adaptint = 10;
   case 'DRAM'
    drscale  = 2;
    adaptint = 10;
  end

  mu   = [0 0];         % center
  cmat = [1 0.9;0.9 1]; % target covariance
  imat = inv(cmat);

  %bpar = [1.3 1.5]; % "bananity" of the target, see bananafun.m
  bpar = [1 1]; % "bananity" of the target, see bananafun.m
  % "sum-of-squares" function as inline function
  bananass = @(x,d) bananafun(x-d.mu,d.bpar,0)*d.imat*bananafun(x-d.mu,d.bpar,0)';
  %bananass = inline('bananafun(x-d.mu,d.bpar,1)*d.imat*bananafun(x-d.mu,d.bpar,1)''','x','d');

  % create input arguments for the dramrun function
  clear model params options

  model.ssfun    = bananass;

  params.par0    = mu; % initial value

  data = struct('mu',mu,'imat',imat,'bpar',bpar);   

  options.nsimu    = nsimu;
  options.adaptint = adaptint;
  options.drscale  = drscale;
  options.qcov     = eye(2)*5; % initial covariance 

  [results,chain] = dramrun(model,data,params,options);
  
  %plot(chain(:,1),chain(:,2),'.','MarkerSize',6)
  
%   figure(1);subplot(2,2,mi);%clf
%   plot(chain(:,1),chain(:,2),'.','MarkerSize',6);
%  
%   %%
%   %% Plot 50% and 95% contours
%   %%
%   c50=1.3863; % critical values from chisq(2) distribution
%   c95=5.9915;
%   hold on
%   [xe,ye]=ellipse(mu,c50*cmat);
%   xyplot(bananafun([xe,ye],bpar,0),'k-','LineWidth',1.5)
%   [xe,ye]=ellipse(mu,c95*cmat);
%   xyplot(bananafun([xe,ye],bpar,0),'k-','LineWidth',1.5)
%   hold off
% 
%   %%% count points inside prob. regions
%   d = mahalanobis(bananafun(chain,bpar,1),mu,imat,1);
%   cc50 = sum(d<c50)/options.nsimu;
%   cc95 = sum(d<c95)/options.nsimu;
% 
%   %%% chain mean and covariance matrix
%   cmean = mean(bananafun(chain,bpar,1));
%   ccov  = cov(bananafun(chain,bpar,1));
% 
%   title(sprintf('%s: %3.1f%% < c50,  %3.1f%% < c95',upper(method),cc50*100,cc95*100))
%   xlabel('\theta_1'); ylabel('\theta_2');
% 
%   axis([-5 5 -5 15]); % fix the axis to comparable between runs
% 
%   % add integrated autocorrelation time
%   text(-1,-8,sprintf('\\tau = %4.1f',mean(iact(chain))));
%   
%   %pngfile(sprintf('banana1%s.png',method));
% 
%   figure(2); clf
%   subplot(2,1,1); plot(chain(:,1),'.')
%   subplot(2,1,2); plot(chain(:,2),'.')
% 
%   %pngfile(sprintf('bananae%s.png',method));
%   
%   figure(3)
%   subplot(4,2,mi)
%   colormap(brewermap([],'Blues'))
%   [bandwidth,density,X,Y] = kde2d(chain,2^7,[-5,-5],[5,15]);
%   imagesc([-5,5],[15,-5],flipud(density))
%   set(gca,'YDir','Normal')
%   hold on
%   bananaplot(0,0,1,1,0.9,1,1,'r--',0.5)
%   bananaplot(0,0,1,1,0.9,1,1,'r--',0.95)
%   text(-4.5,-3,methods{mi})
  
end

passed_chain = chain;
