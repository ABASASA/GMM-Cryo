% preliminaries
fprintf('\n  \n  \n  \n')
clear;
rng(12)
%% generating "ground truth"
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('loading data..')
P = 3;% Bandlimit cutoff
beta  = 1;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy
gridSize = 23;  % number of voxels in each dimenison. take odd.
total_N      = 20000; % Number of projection
numberOfExpriemntsRepet = 64; % number of createing projections' sets (for chainging random dist)
numOfIntialGuesses = 3; % Number of intial guesses
numMaxIter = 100; % Number of Optimizatio iterations

isPadding = true; % to truncate the volume
numOfStepsGMM = 1; % Number of GMM steps

duplicateExpFlag = false; % repreat expirements


total_Ns = total_N * ones(numberOfExpriemntsRepet,1); % the number of projections in each expirements
svPercentArray = [1e-5 * ones(numberOfExpriemntsRepet,1)]; % Optinal for futuere to add parameter to W's computation
sigmas = [0.05* ones(numberOfExpriemntsRepet,1)]; % Noise sigma

% If you want to preform the same tests for 2 sets of params, i.e., same
% projecitons set. Then duplicateExpFlag = True; and enter the relveant
% values to "sigmas" below (or add "total_Ns"):
if duplicateExpFlag
    sigmas = [0.05* ones(numberOfExpriemntsRepet/2,1); 1* ones(numberOfExpriemntsRepet/2,1)];
end
filepathFather = 'ResultsGMM/TestB_seed_12_20000_Ob_3_Guess_sigma005-120521/'; % Path to save results
%% Simulate kspace function

% Sum of Gaussian
sigma = 200;


T     = (gridSize/5)*[0, 0, 0;
                      0.15, 0.10, 0.10;
                      -0.1, -0.15, 0;
                      -0.1, 0.05, -0.05;
                      0, 0.1, -0.1;]';
g     =  @(k) exp(1 * 1i*k'*T(:,1)) .* exp(-0.7   * pi^2/sigma*sum(k.*k)).' + ...
              exp(1 * 1i*k'*T(:,2)) .* exp(-0.5 * pi^2/sigma*sum(k.*k)).' + ...
              exp(1 * 1i*k'*T(:,3)) .* exp(-0.4   * pi^2/sigma*sum(k.*k)).' + ...
              exp(1 * 1i*k'*T(:,4)) .* exp(-0.6   * pi^2/sigma*sum(k.*k)).' + ...
              exp(1 * 1i*k'*T(:,5)) .* exp(-0.8   * pi^2/sigma*sum(k.*k)).' ;


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

% Trncated
if isPadding
    [gamma,A] = gamma_truncate_2(gamma,AFull);
else
    A = AFull;
end

%% calculating moments
fprintf('Calculating moments: \n');

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

fprintf('Run preprocessing...');
tic
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));


%% naming
nameit = ['balls_example_grid_size_',num2str(gridSize)];
tt     = datetime('now');
nameit = [nameit,'P_',num2str(P),'_',num2str(tt.Hour),'_',num2str(tt.Minute)];

%% Create Expiriments
dataStruct = cell(numberOfExpriemntsRepet,1); % multpy by 2 because we have 2 algorithms

%% distribution

load ('SO3_fifteen.mat');

BsCell = cell(numberOfExpriemntsRepet,1);
guessesArray = cell(numberOfExpriemntsRepet,1);
if duplicateExpFlag
    BsCell = cell(numberOfExpriemntsRepet,1);
    guessesArray = cell(numberOfExpriemntsRepet,1);
    for i = 1 : numberOfExpriemntsRepet / 2
       B = get_B_inplane(P);

        [AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
        [BB,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
        scl = 1/B{1};
        for j=1:numel(B)
            B{j} =  scl*B{j};
        end 
        BsCell{i} = B;
        BsCell{i + numberOfExpriemntsRepet/2} = B;
% %         BsCell{i + 2 * numberOfExpriemntsRepet/3} = B;


        intialGuessCells = {};
        vec_AB_GT = A_B_to_vec_inplane(A, B, M-1);
        vec_A_GT  = A_B_to_VecAB(A, [], M-1);
        for j = 1 : numOfIntialGuesses
            if P == 1
                intialGuessCells{j} = -vec_AB_GT + rand(size(vec_AB_GT))*0.001;
            else
                intialGuessCells{j} = get_initial_guess_inplane(P, AI, AE, Acon, size(vec_AB_GT));
            end
        end
        guessesArray{i} = intialGuessCells;
        guessesArray{i + numberOfExpriemntsRepet/2} = intialGuessCells;
        %guessesArray{i + 2 * numberOfExpriemntsRepet/3} = intialGuessCells;
    end
end
 

%% Createt data for optimization
parfor iExp = 1 : numberOfExpriemntsRepet
    total_N = total_Ns(iExp);
    svPercent = svPercentArray(iExp);
    sigmaScalar = sigmas(iExp);
%   
    if  duplicateExpFlag
        B = BsCell{iExp};
        intialGuessCells = guessesArray{iExp}
    else
        B = get_B_inplane(P);

        [AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
        [BB,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
        scl = 1/B{1};
        for i=1:numel(B)
            B{i} =  scl*B{i};
        end

        intialGuessCells = {};
        % vectorized version
        vec_AB_GT = A_B_to_vec_inplane(A, B, M-1);
        vec_A_GT  = A_B_to_VecAB(A, [], M-1);
        for i = 1 : numOfIntialGuesses
            if P == 1
                intialGuessCells{i} = -vec_AB_GT + rand(size(vec_AB_GT))*0.001;
            else
                intialGuessCells{i} = get_initial_guess_inplane(P, AI, AE, Acon, size(vec_AB_GT));
            end
        end
% %     
    end
    
    fprintf('DONE \n');

    
    %% Define noise params (Sigma)
    [sigmaNoise, cov_c] = WhiteNoiseVec(gridSize,floor((gridSize)/2), sigmaScalar);
%     [sigmaNoise, cov_c] = Noise3on3GassianVec(gridSize,floor((gridSize)/2), sigmaScalar);

    %% Genetare observations
    [projs, proj_PSWF, weight, ~, SNR] = GenerateObservationsPSWF(total_N,...
        gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d, sigmaNoise);
    


    %% Compute bias for second moments           
    [Bias] = ComputeBiasSecondMomentAnalyticly(gridSize, floor(gridSize/2), beta,...
                                                eps_p, [], cov_c, gamma);
    %% Estimate m1_hat and m2_hat
    [m1_hat, m2_hat] = EstimateMoments(proj_PSWF, weight, total_N, gamma, Bias);

    %% Create trimming data
    m1Size = size(m1_hat);
    m2Size = size(m2_hat);

    [vecBoolMoments, ~, ~] = TrimMoments_inplane(m1Size, m2Size, gamma, L, P);
    
    
    %% Reg tech
   % [sigmaNoise, cov_c] = WhiteNoiseVec(gridSize,floor((gridSize)/2), 0.05);

    %noise = reshape(mvnrnd(zeros(size(sigmaNoise,1),1), sigmaNoise, total_N)', [gridSize, gridSize,total_N] );
    %projs = projs + noise;
    %[proj_PSWF] = Projection2PSWF(projs, beta, eps_p, gamma);
    %[Bias] = ComputeBiasSecondMomentAnalyticly(gridSize, floor(gridSize/2), beta,...
     %                                           eps_p, [], cov_c, gamma);
    %% Compute W for GMM
%   Ideal
     [W, ~,OmegaCut] = ComputeW_inplane(A, B, proj_PSWF, weight, total_N,...
                                 gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments, svPercent, Bias);
%   ID
%    W = eye(sum(vecBoolMoments)); % Define Id matrix
    fprintf('Done  \n');

    %% setting the optimization
    fprintf('Setting the optimization \n');
    
    N = length(intialGuessCells{1});
    dataStruct{iExp,1}.W = W;
    dataStruct{iExp,1}.intialGuessCells = intialGuessCells;
    dataStruct{iExp,1}.N = N;
    dataStruct{iExp,1}.m1_hat = m1_hat;
    dataStruct{iExp,1}.m2_hat = m2_hat;
    dataStruct{iExp,1}.vecBoolMoments = vecBoolMoments;
    dataStruct{iExp,1}.weight = weight;
    dataStruct{iExp,1}.svPercent = svPercent;
    dataStruct{iExp,1}.SNR  = SNR;
    dataStruct{iExp,1}.Bias = Bias;
    dataStruct{iExp,1}.proj_PSWF = proj_PSWF;
    dataStruct{iExp,1}.B = B;
end
%% running optimization

ansCell = cell(numberOfExpriemntsRepet *  2 * numOfIntialGuesses, 1); % cell  to save results
% delete(gcp('nocreate'))

parfor iOpt = 1 : 2  * numOfIntialGuesses * numberOfExpriemntsRepet
    [iOptRow, iOprCol] = ind2sub([numberOfExpriemntsRepet, 2 * numOfIntialGuesses], iOpt);
    W = dataStruct{iOptRow,1}.W;
    N = dataStruct{iOptRow,1}.N;
    m1_hat = dataStruct{iOptRow,1}.m1_hat;
    m2_hat = dataStruct{iOptRow,1}.m2_hat;
    SNR = dataStruct{iOptRow,1}.SNR;
    vecBoolMoments = dataStruct{iOptRow,1}.vecBoolMoments;
    Bias =  dataStruct{iOptRow,1}.Bias;
    intialGuessCells = dataStruct{iOptRow,1}.intialGuessCells;
    svPercent = dataStruct{iOptRow,1}.svPercent;
    proj_PSWF = dataStruct{iOptRow,1}.proj_PSWF;
%     proj_PSWF = [];

    weight = dataStruct{iOptRow,1}.weight;

    if iOprCol > numOfIntialGuesses % GMM

        fprintf(['Run GMM RowIndex: ', num2str(iOptRow) ' Colindex:' num2str(iOprCol - numOfIntialGuesses) '...  \n']);

        initial_guess = intialGuessCells{iOprCol - numOfIntialGuesses};
        %% Run GMM optimization
        [t_optGMM, xEstGMM, xcostGMM, infoGMM] = Iterative_GMM_Inplane(N, P, gamma, C_tensor,...
                        Gamma_mat, sign_mat, m1_hat, m2_hat, vecBoolMoments, proj_PSWF, weight,...
                        svPercent, initial_guess, numMaxIter, numOfStepsGMM, W, Bias);
        fprintf(['Done GMM RowIndex: ', num2str(iOptRow) ' Colindex:' num2str(iOprCol - numOfIntialGuesses)]);
        % Save data GMM

        stGMM = struct;
        stGMM.t_optGMM = t_optGMM;
        stGMM.xEstGMM = xEstGMM;
        stGMM.xcostGMM = xcostGMM;
        stGMM.infoGMM = infoGMM;
        ansCell{iOpt} = stGMM;
    else %LS
        fprintf(['Run LS RowIndex: ', num2str(iOptRow) ' Colindex:' num2str(iOprCol) '...  \n']);
        initial_guess = intialGuessCells{iOprCol};
        %% Run LS optimization
        weightsLs = norm(m1_hat(:))^2/norm(m2_hat(:))^2;
        [t_optLS, xLS, xcostLS, infoLS ] = Opimization_LS_inplane(N, weightsLs, P, gamma, C_tensor, ...            
                Gamma_mat, sign_mat, m1_hat, m2_hat, initial_guess, numMaxIter);
                       %  svPercent, initial_guess, numMaxIter, numOfStepsGMM, W, Bias);
        % Save data LS
        stLS = struct;
        stLS.t_optLS = t_optLS;
        stLS.xLS = xLS;
        stLS.xcostLS = xcostLS;
        stLS.infoLS = infoLS;
        ansCell{iOpt} = stLS;
        fprintf(['Done LS index:' num2str(iOpt)]);



    end

    ansCell{iOpt}.SNR = SNR;

end
fprintf('Finished...  \n');
ansCellMatrix = reshape(ansCell,[numberOfExpriemntsRepet, 2 * numOfIntialGuesses]);

%%
for iExp = 1 : numberOfExpriemntsRepet
    dataStruct{iExp,1}.proj_PSWF =  [];
    dataStruct{iExp,1}.Bias =  [];
    dataStruct{iExp,1}.weight =  [];
end
%% save
mkdir(filepathFather);
save([filepathFather, 'data.mat'], 'ansCellMatrix', 'filepathFather', 'total_Ns', 'beta','total_N', 'numberOfExpriemntsRepet'...
            ,'eps_p','delta', 'svPercentArray', 'isPadding', 'sigmas', 'dataStruct', 'gamma');
% delete(gcp('nocreate'))
% parpool();
%% Analyize results 
for iExp = 1 : numberOfExpriemntsRepet
minIndexLS = inf;
currentCostLS = inf;
minIndexGMM = inf;
currentCostGMM = inf;
% Find intial guess which its optimization's output scored the lowest.
for iOpt = 1: numOfIntialGuesses
    if currentCostLS > ansCellMatrix{iExp, iOpt}.xcostLS
        minIndexLS = iOpt;
        currentCostLS = ansCellMatrix{iExp, iOpt}.xcostLS;
    end
    
    if currentCostGMM > ansCellMatrix{iExp, iOpt + numOfIntialGuesses}.xcostGMM
        minIndexGMM = iOpt + numOfIntialGuesses;
        currentCostGMM = ansCellMatrix{iExp, iOpt + numOfIntialGuesses}.xcostGMM;

    end
    
end
ansCellLS = ansCellMatrix{iExp, minIndexLS};
ansCellGMM = ansCellMatrix{iExp, minIndexGMM};

[A_LS, ~] = VecAB_inplane_to_A_B(ansCellLS.xLS, gamma);
[A_GMM, ~] = VecAB_inplane_to_A_B(ansCellGMM.xEstGMM, gamma);
 
SNR = ansCellMatrix{iExp,1}.SNR;

% if ~isPadding
%     AFull = []; % if is padding false
% end
disp('Change ME')
filepath = [filepathFather, 'Exp' num2str(iExp) '/'];
mkdir(filepath);
%% LS - Align, Compute FSC and rellative error
estimationMethodNameString = 'LS';

[rel_errLS, resAEstFSCLS, figFSCLS] = AnalytsisResultsFCS_Alignment_relError(A, A_LS, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit,c);
rel_errLS        


figFSCGMM = figure;
%% GMM - Align, Compute FSC and rellative error
% subplot(2,1,1);
estimationMethodNameString = 'GMM';
[rel_errGMM, resAEstFSCGMM, figFSCGMM] = AnalytsisResultsFCS_Alignment_relError(A, A_GMM, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit,c, figFSCGMM);
rel_errGMM

save([filepath, 'data.mat'],'rel_errLS','resAEstFSCLS','figFSCLS','rel_errGMM','resAEstFSCGMM',...
        'figFSCGMM', 'SNR');
disp(['LS (FSC, Error) : (' num2str(resAEstFSCLS) ', ' num2str(rel_errLS) ') GMM: (' num2str(resAEstFSCGMM)...
    ', ' num2str(rel_errGMM) ')'])
     
end