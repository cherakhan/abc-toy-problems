error = (0:0.1:10)';

e1 = 0.3;
e2 = 3;

sisson = zeros(length(error),1);
mahala = zeros(length(error),1);

for ii = 1:length(error)
    
   sisson(ii) = exp(-0.5*(error(ii)/e2)*1*(error(ii)/e2));
   
   mahala(ii) = exp(-0.5*(error(ii))*1/e2*(error(ii)));
    
end

plot(error,sisson,error,mahala)