function [gcurrent] = CreateElipsoidVol(centersMatrix, sizeFactorArray)

N = size(centersMatrix, 2);
gs = cell(N,1);
gcurrent = @(k) 0;

for i = 1 : N
    factor = sizeFactorArray(i);
    tmpG = @(k) exp(1i*k'*centersMatrix(:,i)).*exp(-factor*pi^2/sigma*sum(k.*k)).';
%     gs{i}
    oldg = @(k) gcurrent(k);
    gcurrent = @(k) oldg(k) + tmpG(k);
end
end