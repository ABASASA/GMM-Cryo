fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth" 
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('Generating data..')

% parameters
P        = 3;   % distribution expansion length
gridSize = 29;   % volume resolution
delta    = .99;
beta     = .9; % between 0 and 1, controls the expansion length
total_N      = 10;

% getting the map -- volume expansion coefficients
if and(gridSize == 65, beta == .6)
    load('A_coefs_beta_6_grid_65.mat');
    eps_p  = 5e-4;    % Prescribed accuracy in images side
    radius = floor(gridSize/2);
    c      = beta*pi*radius; 
else
    map    = ReadMRC('emd_0409_411.mrc'); % upload the 80X80X80 vol
    vol    = cryo_downsample(map,[gridSize gridSize gridSize]); % Downsamples map
    eps_p  = 5e-4;    % Prescribed accuracy in images side
    radius = floor(gridSize/2);
    c      = beta*pi*radius;              % Nyquist bandlimit, if beta == 1
    % volume and change of coordinates coefficients
    A      = pswf_t_f_3d(vol, beta, delta);
    gamma  = PSWF_2D_3D_T_mat(c, delta, eps_p);
%     save('A_coefs_beta','A','gamma');
end
 
%% distribution
B = get_B_inplane(P);


%% A
fprintf('Calculating 3D coefficients...');
AFull = pswf_t_f_3d(vol, beta, delta);
volBack      = pswf_t_b_3d(AFull,     gridSize, beta, delta);
VisualVol(vol, ['vol_3NirAsaf']);
VisualVol(volBack, ['vol_4NirAsaf']);
fprintf('DONE \n');

%% pre-calculating moments
fprintf('Calculating moments: \n');

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
% P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

fprintf('Run preprocessing...');
tic
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

fprintf('DONE in about %d seconds \n', round(tt));

% first moment
[ m1_true ] = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
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
Ns = floor(logspace(1,2,2));

repNumber = 10;

% end

errorM1 = zeros(repNumber, length(Ns));
errorM2 = zeros(repNumber, length(Ns));

[~, proj_PSWF, weight] = GenerateObservationsPSWFFromVol(volBack,...
                                total_N, gridSize, B, beta, eps_p, gamma);
%%
total_N= Ns(end);

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