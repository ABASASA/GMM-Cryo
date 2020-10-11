    % preliminaries
    fprintf('\n  \n  \n  \n')
    clear;

    %% generating "ground truth"
    fprintf('**** Starting simulation ****\n');
    fprintf('\n  \n')
    fprintf('loading data..')
    P  = 3;     % aribitrary, chosen to be small
    B = get_B_inplane(P); % Random inplane dist
    beta  = 1;       % Bandlimit ratio
    delta = 0.99;    % Truncation parameter
    eps_p = 1e-3;    % Prescribed accuracy
    totalN = 10000;%1e5;
    repNumber = 5; % number of repeats
    paramsName = 'C1_params';

    %% Simulate kspace function
    gridSize = 15;  % number of voxels in each dimenison. take odd.
    sigmaMoment = WhiteNoise(gridSize, 0);         % Enter 0 to work without noise.

%% uncomment for my vol
%     % Sum of Gaussian
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
%     paramsName = 'C1_params';
%     
%     vol    = cryo_gaussian_phantom_3d(paramsName,gridSize,1);
%     volvec = vol(:);
%     treshold = 0;
%     volvec(volvec < treshold) = 0;
%     vol = reshape(volvec, size(vol));
%     vol  = vol*floor(size(x_2d,1)/2);   % Nir Sunday ?

    %% calculate 3D expansion and gamma coefficients


    fprintf('Calculating 3D coefficients...');
    AFull = pswf_t_f_3d(vol, beta, delta);
    fprintf('DONE \n');

    fprintf('Calculating gamma coefficients...');
    radius = floor(gridSize/2);
    c     = beta*pi*radius;              % Nyquist bandlimit
    gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);
    %Trncated
%     [gamma,A] = gamma_truncate_2(gamma,AFull);
    A = AFull;
    fprintf('DONE \n');
%     vol = pswf_t_b_3d(A,  gridSize, c, delta);% *floor(size(x_2d,1)/2);
%     A = pswf_t_f_3d(vol, beta, delta);

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
%     [ m1_true ] = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);

    fprintf('First moment is DONE. \n');

    % second moment
    fprintf('Mu2 calculation...');
    tic
    [m2_true] = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
    tt = toc;
    fprintf('DONE in about %d seconds \n', round(tt));

    %% Compute Projecitons
%     paramsName = [];%'C1_params' ;
    % Ns = floor(linspace(5,10000,50));
    Ns = floor(logspace(1,log10(totalN),20));


    % end

    errorM1 = zeros(repNumber, length(Ns));
    errorM2 = zeros(repNumber, length(Ns));

    total_N= Ns(end);
    proj_Cell = struct;

    %% Create projeciton for the number of reptitions (in order to parallel)
    for iRep = 1 : repNumber
        currentstruct = struct;
        tic
%         [~, proj_PSWF, weight, SNR] = GenerateObservationsPSWF(paramsName,...
%                                     Ns(end), gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d, sigmaMoment);
        [~, proj_PSWF, weight, SNR] =  GenerateObservationsPSWFFromAspire(paramsName,...
                                    Ns(end), gridSize, B, beta, eps_p, gamma, sigmaMoment, vol);
        proj_Cell(iRep).proj_PSWF = proj_PSWF;
        proj_Cell(iRep).weight = weight;
        proj_Cell(iRep).sigmaMoment = sigmaMoment;
        toc
    end
    % Compute bias
    [Bias] = ComputeBiasSecondMomentAnalyticly(gridSize,  floor(gridSize/2), beta,...
                                                eps_p, [], sigmaMoment, gamma);
    Bias = 0;
    %% compute the empricial moments
    for iRep = 1 : repNumber
        total_N = Ns(end);
        tmpErrorM1 = zeros(length(Ns),1);
        tmpErrorM2 = zeros(length(Ns),1);
        fprintf('This is the %d round \n', iRep);
    %         perm_list = randperm(size(R,3));
        perm_list = randperm(total_N);

        proj_PSWF = proj_Cell(iRep).proj_PSWF;
        weight = proj_Cell(iRep).weight;
        % Add ratio to the sigma:

        parfor indexN = 1 : length(Ns)
            currentPer = perm_list(1:Ns(indexN));
            %% Estimate m1_hat and m2_hat
            [m1_hat, m2_hat] = EstimateMoments(proj_PSWF(:,:,currentPer), weight(currentPer), Ns(indexN), gamma, Bias);
            tmpErrorM1(indexN) = norm(m1_hat(:) - m1_true(:)).^2 ./ (norm(m1_true(:))^2);
            tmpErrorM2(indexN) = norm(m2_hat(:) - m2_true(:)).^2 ./ (norm(m2_true(:))^2);
        end
        errorM1(iRep,:) = tmpErrorM1;
        errorM2(iRep,:) = tmpErrorM2;

    end
    errorM1 = errorM1.';
    errorM2 = errorM2.';


    %% Display
    fig = figure;
    subplot(2,1,1);
    loglog(Ns, mean(errorM1,2),'*-');
    title('Convargance for M1');

    subplot(2,1,2);
    loglog(Ns, mean(errorM2,2),'*-');
    title('Convargance for M2');
    %%
    % filepath = 'ResultsGMM/020920- ConvargeOfEstimator/Guassains-P_2_Beta_1_WithLine/';
    % flagInput = input(['To Save in: ',filepath ]);
    % if flagInput
    %     mkdir(filepath);
    %     saveas(fig, [filepath,'ConvargeFig'],'fig');
    %     saveas(fig, [filepath,'ConvargeFig'],'jpg');
    %     save([filepath,'data.mat'], 'beta','P','errorM1', 'errorM2', 'vol');
    % end
