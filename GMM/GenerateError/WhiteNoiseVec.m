function [covaraince, cov_c] = WhiteNoiseVec(gridSize,L, sigmaScalar)

%% Generate regular grid and basis functions (if needed)
x_1d_grid = -L:1:L;   % - Odd number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= L);


sigma = sigmaScalar * ones(gridSize, gridSize);


%% Cut

covaraince = sigma(:).^2;
cov_c = diag(covaraince(points_inside_the_circle));
covaraince(~points_inside_the_circle) = 0;
covaraince = diag(covaraince);

end