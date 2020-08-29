function [projs, proj_PSWF, weight] = GenerateObservationsPSWFFromVol(vol,...
                                total_N, gridSize, B, beta, eps_p, gamma)

[projs ,weight] = ComputeProjectionFromVol(total_N, gridSize, vol, B);

[proj_PSWF] = Projection2PSWF(projs, beta, eps_p, gamma);
end