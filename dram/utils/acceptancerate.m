function [AR] = acceptancerate(chain)
% Calculate acceptance rate for 5 par chain

counter = 0;

for ii = 2:length(chain)
    if chain(ii-1,1)~=chain(ii,1) || chain(ii-1,2)~=chain(ii,2) || chain(ii-1,3)~=chain(ii,3) || chain(ii-1,4)~=chain(ii,4) || chain(ii-1,5)~=chain(ii,5)
        counter = counter+1;
    end
end

AR = counter/length(chain);