% preliminaries
fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth"
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('loading data..')
P  = 2;     % aribitrary, chosen to be small
B = get_B_inplane(P);



%% Simulate kspace function
gridSize = 19;  % number of voxels in each dimenison. take odd.

% Sum of Gaussian
% vol1    = cryo_gaussian_phantom_3d('C1_params',gridSize,1);
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
beta  = 1;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy

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

fprintf('Moments calculation...');
% first moment
[ m1_true ] = FirstMoment_PSWF_v1_fixed(A, B, gamma);%Gamma_mat, sign_mat);
fprintf('First moment is DONE. \n');

% second moment
fprintf('Mu2 calculation...');
tic
[m2_true] = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

%% Compute Projecitons
paramsName = [];%'C1_params' ;
Ns = floor(linspace(5,10000,50));
Ns = floor(logspace(1,5,20));

repNumber = 10;

% end

errorM1 = zeros(repNumber, length(Ns));
errorM2 = zeros(repNumber, length(Ns));
[~, proj_PSWF, weight] = GenerateObservationsPSWF(paramsName,...
                                Ns(end), gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d);
total_N= Ns(end);

%%
for iRep = 1 : repNumber
    total_N = Ns(end);
    tmpErrorM1 = zeros(length(Ns),1);
    tmpErrorM2 = zeros(length(Ns),1);
    fprintf('This is the %d round \n', iRep);
%         perm_list = randperm(size(R,3));
    perm_list = randperm(total_N);

    
    parfor indexN = 1 : length(Ns)
        currentPer = perm_list(1:Ns(indexN));
        %% Estimate m1_hat and m2_hat
        [m1_hat, m2_hat] = EstimateMoments(proj_PSWF(:,:,currentPer), weight(currentPer), Ns(indexN), gamma);
        tmpErrorM1(indexN) = norm(m1_hat(:) - m1_true(:)).^2 ./ (norm(m1_true(:))^2);
        tmpErrorM2(indexN) = norm(m2_hat(:) - m2_true(:)).^2 ./ (norm(m2_true(:))^2);
    end
    errorM1(iRep,:) = tmpErrorM1;
    errorM2(iRep,:) = tmpErrorM2;

end
errorM1 = errorM1.';
errorM2 = errorM2.';


%% Display
figure;
subplot(2,1,1);
loglog(Ns, mean(errorM1,2),'*-');
title('Convargance for M1');

subplot(2,1,2);
loglog(Ns, mean(errorM2,2),'*-');
title('Convargance for M2');
