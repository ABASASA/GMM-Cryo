function [m2Bias] = ComputeBiasSecondMomentAnalyticly(sizeImages, L, beta,...
                                                eps_p, PSWF_Nn_p, sigma, gamma)
%% Generate regular grid and basis functions (if needed)
x_1d_grid = -L:1:L;   % - Odd number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= L);

if isempty(PSWF_Nn_p)
    x = x_2d_grid(points_inside_the_circle);
    y = y_2d_grid(points_inside_the_circle);
    [ PSWF_Nn_p.samples, PSWF_Nn_p.alpha_Nn, PSWF_Nn_p.ang_freq, PSWF_Nn_p.rad_freq ] = PSWF_gen_v3( L+1, L, beta, eps, eps_p, x, y );
    PSWF_Nn_p.L = L;
    PSWF_Nn_p.beta = beta;
    PSWF_Nn_p.T = eps_p;
end

%% compute expansion coefficients
nImages = 1;
images_c = ones(sizeImages * sizeImages, 1);    % v
images_c = images_c(points_inside_the_circle(:),:);     % Take points inside the unit disk

sigma = sigma(:);
sigma_c = sigma(points_inside_the_circle(:));

varSigma = diag(sigma_c .^ 2);
% Idsigma = sigma^2 * eye(length(images_c), length(images_c));

PSWF_m = bsxfun(@times,PSWF_Nn_p.samples,abs((beta/2)*PSWF_Nn_p.alpha_Nn.').^2); 
PSWF_mTans = PSWF_m';

    
m2Bias  = zeros(max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0), ...
                max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0));
for m1=0:max(gamma.ang_idx_2d)
    for k1=0:nnz(gamma.ang_idx_2d==m1)-1
         is = PSWF_Nn_p.ang_freq== (m1);
         i = find(cumsum(is) == (k1 + 1),1);
        for m2=0:max(gamma.ang_idx_2d)  
            for k2=0:nnz(gamma.ang_idx_2d==m2)-1
                js = PSWF_Nn_p.ang_freq== (m2);
                j = find(cumsum(js) == (k2 + 1),1);
                
                

                m2Bias(m1 + 1, k1 + 1, m2 + 1, k2 + 1) = PSWF_mTans(i,:) * varSigma * PSWF_mTans(j,:)';
            end
        end
    end
end
end