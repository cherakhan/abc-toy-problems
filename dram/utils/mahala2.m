function d=mahala2(x,mu,CM)
% x and mu are row vectors

% command that worked: md = sqrt((x-mu)*Isigma*(x-mu)')
%  md = sqrt(([0,1]-mu)*Isigma*([0,1]-mu)')

iCM = inv(CM+(eye(length(CM))*10e-10));

d = sqrt((x-mu)*iCM*(x-mu)');