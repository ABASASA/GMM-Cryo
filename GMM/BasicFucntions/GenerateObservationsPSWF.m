function [projs, proj_PSWF, weight, sigmaMat, SNR] = GenerateObservationsPSWF(paramsName,...
                                total_N, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d, sigmaMat)

[projs ,weight] = ComputeProjection (total_N, gridSize, paramsName, B, g, x_2d, y_2d);


if sum(sigmaMat) ~=0
%     [projs, ~, ~, sigma] = cryo_addnoise(projs, SNR,'gaussian');
    
    SNR = mean(vec(sum(projs,3) / total_N)) / mean(sigmaMat(:));

    
    projs = projs + randn(size(projs)) .* sigmaMat;
else
    SNR = 0;
end
[proj_PSWF] = Projection2PSWF(projs, beta, eps_p, gamma);


end