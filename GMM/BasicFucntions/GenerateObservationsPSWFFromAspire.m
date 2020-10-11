function [projs, proj_PSWF, weight, sigmaMat, SNR] = GenerateObservationsPSWFFromAspire(paramsName,...
                                total_N, gridSize, B, beta, eps_p, gamma, sigmaMat, vol)

[projs ,weight] = ComputeProjectionFromAspire (vol, total_N, gridSize, paramsName, B);


if sum(sigmaMat) ~=0
%     [projs, ~, ~, sigma] = cryo_addnoise(projs, SNR,'gaussian');
    
    SNR = mean(vec(sum(projs,3) / total_N)) / mean(sigmaMat(:));

    
    projs = projs + randn(size(projs)) .* sigmaMat;
else
    SNR = 0;
end
[proj_PSWF] = Projection2PSWF(projs, beta, eps_p, gamma);


end