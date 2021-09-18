function [projs, proj_PSWF, weight, sigmaMat, SNR] = GenerateObservationsPSWF(...
                                total_N, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d, sigmaMat)

[projs ,weight] = ComputeProjection (total_N, B, g, x_2d, y_2d);


if sum(abs(sigmaMat(:))) ~=0
    
    noise = reshape(mvnrnd(zeros(size(sigmaMat,1),1), sigmaMat, total_N)', [gridSize, gridSize,total_N] );
    SNR = mean(vec(sum(projs.^2,3)/ total_N)) / mean(vec(sum(noise.^2,3)/ total_N));
    

    projs = projs + noise;
    
else
    SNR = 0;
end
[proj_PSWF] = Projection2PSWF(projs, beta, eps_p, gamma);


end