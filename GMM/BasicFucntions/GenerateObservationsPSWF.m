function [projs, proj_PSWF, weight] = GenerateObservationsPSWF(paramsName,...
                                total_N, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d)

[projs ,weight] = ComputeProjection (total_N, gridSize, paramsName, B, g, x_2d, y_2d);

[proj_PSWF] = Projection2PSWF(projs, beta, eps_p, gamma);
end