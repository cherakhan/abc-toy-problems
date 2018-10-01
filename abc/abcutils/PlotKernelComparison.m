function [] = PlotKernelComparison()

colors = brewermap(4,'Set1');

error = 0:0.01:10;
tolerance = 1;

g = zeros(length(error),1);
u = zeros(length(error),1);

for ii = 1:length(error)
    g(ii) = exp(-0.5*gaussiankernel(error(ii),tolerance));
    u(ii) = exp(-0.5*uniformkernel(error(ii),tolerance));
end
plot(error,g,'color',colors(2,:))
hold on
plot(error,u,'color',colors(4,:))
legend({'Gaussian kernel','Uniform kernel'})
legend boxoff

end