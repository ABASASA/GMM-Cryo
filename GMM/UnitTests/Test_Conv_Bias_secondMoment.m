% preliminaries
fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth"
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('loading data..')
P  = 3;     % aribitrary, chosen to be small
B = get_B_inplane(P);
beta  = 1;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy

%% Simulate kspace function
gridSize = 19;  % number of voxels in each dimenison. take odd.

sigma = 200;

T     = (gridSize/5)*[0 0 0; 0.08 0.1 0; -0.1 0 0.1]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
g     =  @(k) exp(1i*k'*T(:,1)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,2)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,3)).*exp(-pi^2/sigma*sum(k.*k)).' ;%+ ...
%% Vol on kspace over cartesian grid

radius   = floor(gridSize/2);
if mod(gridSize,2)==0
    x_1d = (-radius:1:radius-1);%/radius;   % - Even number of points
else
    x_1d = (-radius:1:radius);%/radius;   % - Odd number of points
end

% the 3d grid
[x_3d,y_3d,z_3d] = meshgrid(x_1d,x_1d,x_1d);
vec_3d_grid      = [x_3d(:),y_3d(:),z_3d(:)].'*pi/2;

% parameterization of the 2d slices
x_2d = x_3d(:,:,1);
y_2d = y_3d(:,:,1);

% evaluate the function
volf = g(vec_3d_grid);

% back to real domain
volk = reshape(volf, gridSize, gridSize, gridSize);
vol  = real(fftshift(ifftn(ifftshift(volk)))); % in real domain
vol  = vol*floor(size(x_2d,1)/2);   % Nir Sunday ?

%% calculate 3D expansion and gamma coefficients


fprintf('Calculating 3D coefficients...');
AFull = pswf_t_f_3d(vol, beta, delta);
fprintf('DONE \n');

fprintf('Calculating gamma coefficients...');
radius = floor(gridSize/2);
c     = beta*pi*radius;              % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);
%Trncated
[gamma,A] = gamma_truncate_2(gamma,AFull);
% A = AFull;
fprintf('DONE \n');

%% calculating moments
fprintf('Calculating moments: \n');

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

fprintf('Run preprocessing...');
tic
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

%% Compute Projecitons
paramsName = [];%'C1_params' ;

sigma = 2.5;
Ns = floor(logspace(0.5,5.5,20));
repNumBias = 5;
Biases = {};

% analBias = ComputeBiasMatrixEmpricly(B, sigma, gridSize,  g, x_2d, y_2d, gamma, 1, repNumBias);
parfor ii = 1 : length(Ns)
 Biases{ii} = ComputeBiasMatrixEmpricly(B, sigma, gridSize,  g, x_2d, y_2d, gamma, Ns(ii), repNumBias);
end
%%
% 
% [~, proj_PSWF, ~, ~] = GenerateObservationsPSWF(paramsName,...
%                                     1, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d, 0);
% % m1Size = size(proj_PSWF);
% m2Size = size(Biases{end});
% m1Size = m2Size(1:2);
% [vecBoolMoments, boolM1, boolM2] = TrimMoments_inplane(m1Size,...
%                         m2Size, gamma, L, P);
% analBias = zeros(size(Biases{1}));
% for m=0:max(gamma.ang_idx_2d)
%     for k=0:nnz(gamma.ang_idx_2d==0)-1
% %         Grid 15
% %         if (boolM1(m+1,k+1) ~= 0)
% %              analBias(m+1,k+1, m +1, k+1) = sigma^2 * 0.02040047;
% %         end
% %         if m == 0 & k == (nnz(gamma.ang_idx_2d==0) - 1)
% %             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 * 0.0128724;
% %         elseif m == max(gamma.ang_idx_2d) & k == (nnz(gamma.ang_idx_2d==0) - 2)
% %             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 *0.0134306;
% %         elseif m == max(gamma.ang_idx_2d)-1 & k == (nnz(gamma.ang_idx_2d==0) - 2)
% %             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 * 0.018947370;
% %         end
% 
% %         Grid 19
%         if (boolM1(m+1,k+1) ~= 0)
%              analBias(m+1,k+1, m +1, k+1) = sigma^2 * 0.0122992099;
%         end
% 
%         if m == 0 & k == (nnz(gamma.ang_idx_2d==0) - 1)
%             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 * 0.008039786982;
% %         elseif m == 0 & k == (nnz(gamma.ang_idx_2d==0) - 2)
% %             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 * 0.008039786982;
%         elseif m == max(gamma.ang_idx_2d) & k == (nnz(gamma.ang_idx_2d==0) - 2)
% %             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 * 0.0134306;
%         elseif m == max(gamma.ang_idx_2d)-1 & k == (nnz(gamma.ang_idx_2d==0) - 2)
%             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 * 0.00826;
%         elseif m == 1 & k == (nnz(gamma.ang_idx_2d==0) - 2)
%             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 * 0.011425929461407;
%         elseif m == max(gamma.ang_idx_2d) & k == (nnz(gamma.ang_idx_2d==0) - 3)
%             analBias(m + 1,k + 1,m + 1,k +1) = sigma^2 * 0.011534171917999;
%         end
%     end
% end

%%
beta  = 1;       % Bandlimit ratio
eps_p = 1e-3;    % Prescribed accuracy
% LL = floor(size(gridSize,1)/2);
[proj ,weight] = ComputeProjection (1, gridSize, paramsName, B, g, x_2d, y_2d);
[analBias] = ComputeBiasSecondMomentAnalyticly(size(proj,1), floor(size(proj,1)/2), beta,...
                                                eps_p, [], sigma, gamma);
                                            
                                            
                                            
1
%%
diff1 = vec(Biases{end} - analBias);
figure;
plot(1: length(diff1),diff1, '*--');


relError = zeros(length(Ns),1);
for jj = 1: length(Ns)
    relError(jj) = norm(vec(Biases{jj} - analBias) , 'fro') / norm(vec(analBias) , 'fro');
end
figure
loglog(Ns(1:end), relError, 'r*-');
xlabel('#projection')
ylabel('Rel.Error')
title('Bias convargance')
