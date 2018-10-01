function [] = bananaplot(muX,muY,varX,varY,correlation,b1,b2,col,interval)

addpath('../dram/utils')
addpath('../dram')

c95 = chiqf_m(interval,2);


mu = [muX,muY];
rho = correlation;
CM = [varX,rho*sqrt(varX*varY);rho*sqrt(varX*varY),varY];

[xe1,ye1] = ellipse(mu,CM*c95);
xyplot(bananafun([xe1,ye1],[b1,b2],1),col,'LineWidth',1.5)